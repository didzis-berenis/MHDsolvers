/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "incompressibleConductingVoFEtransfer.H"
#include "fvcGrad.H"
#include "fvcDiv.H"
#include "fvcSurfaceIntegrate.H"
#include "fvmLaplacian.H"
#include "surfaceInterpolate.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include <map>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(incompressibleConductingVoFEtransfer, 0);
    addToRunTimeSelectionTable
    (
        solver,
        incompressibleConductingVoFEtransfer,
        fvMesh
    );
}
}


// * * * * * * * * * * * * * * Electric-field transfer * * * * * * * * * * * //
//
// Optional coupling mode, selected in constant/electromagneticModel:
//
//     receivedField       J;     // J (default) | E
//     eTransferProjection none;  // none (default) | full | incremental
//
// receivedField J (or the key absent): the received fields are the current
// density, and the code below is not used. The J-transfer path is unchanged.
//
// receivedField E: the fields received into Jre/Jim at a coupling exchange
// are the harmonic electric field (real and imaginary part), not the current
// density. The module then forms the current itself from the CURRENT alpha1
// and U at every step, with E and B held between exchanges:
//
//     W      = E + (U - U_exchange) x B
//     sigma  = a sigma_c,  sigmaS = (floor + (1 - floor) a) sigma_c,
//     a      = clip(alpha1, 0, 1)
//
// The exported E already contains the motional term of the exchange,
// U_exchange x B (the field solver adds v x curl(A) to E when a Lorentz
// velocity is given), so only the velocity change since the exchange is
// added. U_exchange is the velocity stored at the exchange (storeU).
//
// eTransferProjection
//   none (default): J = sigma W, no potential correction.
//   full:           face flux F = (sigmaS W)_f . S_f - flux(phi), with phi from
//                   laplacian(sigmaS, phi) = div((sigmaS W)_f . S_f), so that
//                   sum_f F = 0 in every cell; F = 0 on non-coupled boundary
//                   patches (zeroGradient phi, source zeroed there);
//                   J = (sigma/sigmaS) (1/V) sum_f F_f (x_f - x_P)
//                   (conservative reconstruction).
//   incremental:    phi from laplacian(sigmaS, phi) = div(psi_now)
//                   - div(psi_exchange); J = sigma W
//                   + (sigma/sigmaS) (1/V) sum_f (-flux(phi))_f (x_f - x_P),
//                   i.e. sigma E at the exchange and only the divergence
//                   created since the exchange projected out.
//
// In every mode:
//   Joule heat JJsigma = 0.5 sigma (|W_re - grad phi_re|^2
//                                 + |W_im - grad phi_im|^2)
//   force      JxB     = 0.5 (J_re x B_re + J_im x B_im)
// phi is copied to PotEre/PotEim and the received E is written as
// eTransferERe/eTransferEIm for diagnostics. The full and incremental modes
// solve for eTransferPhiRe/eTransferPhiIm, which need an fvSolution entry.
//
// The state is kept at file scope (per mesh) and in the mesh registry, so
// that this module adds no data members to incompressibleConductingVoF.

namespace
{
    //- Lower bound of the solve conductivity, (floor + (1 - floor) a) sigma;
    //  mirrors the private constant floor_ of incompressibleConductingVoF
    const Foam::scalar eTransferFloor = 1e-6;

    //- Received-field bookkeeping of the E path, per mesh
    struct eTransferState
    {
        Foam::label version = 0;        // incremented at every received field
        Foam::label doneVersion = -1;   // version of the last projection
        Foam::label doneTimeIndex = -1; // time index of the last projection
        Foam::label refVersion = -1;    // version of the incremental reference
    };

    eTransferState& eState(const Foam::fvMesh& mesh)
    {
        static std::map<Foam::word, eTransferState> states;
        return states[mesh.name()];
    }

    Foam::volScalarField& eRegScalar
    (
        const Foam::fvMesh& mesh,
        const Foam::word& name,
        const Foam::dimensionSet& dims,
        const bool write
    )
    {
        using namespace Foam;
        if (!mesh.foundObject<volScalarField>(name))
        {
            volScalarField* fPtr
            (
                new volScalarField
                (
                    IOobject
                    (
                        name,
                        mesh.time().name(),
                        mesh,
                        IOobject::NO_READ,
                        write ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensionedScalar(dims, 0),
                    zeroGradientFvPatchScalarField::typeName
                )
            );
            fPtr->store(fPtr);
        }
        return mesh.lookupObjectRef<volScalarField>(name);
    }

    Foam::volVectorField& eRegVector
    (
        const Foam::fvMesh& mesh,
        const Foam::word& name,
        const Foam::dimensionSet& dims,
        const bool write
    )
    {
        using namespace Foam;
        if (!mesh.foundObject<volVectorField>(name))
        {
            volVectorField* fPtr
            (
                new volVectorField
                (
                    IOobject
                    (
                        name,
                        mesh.time().name(),
                        mesh,
                        IOobject::NO_READ,
                        write ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensionedVector(dims, Zero),
                    zeroGradientFvPatchVectorField::typeName
                )
            );
            fPtr->store(fPtr);
        }
        return mesh.lookupObjectRef<volVectorField>(name);
    }

    //- JxB and JJsigma of the E path from the current J, B and the stored heat
    void eTransferForceAndHeat
    (
        const Foam::fvMesh& mesh,
        Foam::electromagneticModel& em
    )
    {
        using namespace Foam;
        volVectorField& JxB = em.getVectorFieldRef("JxB");
        JxB.primitiveFieldRef() =
            0.5
           *(
                (em.J().primitiveField() ^ em.B().primitiveField())
              + (em.J(true).primitiveField() ^ em.B(true).primitiveField())
            );
        JxB.correctBoundaryConditions();
        volScalarField& JJsigma = em.getScalarFieldRef("JJsigma");
        JJsigma.primitiveFieldRef() =
            mesh.lookupObject<volScalarField>("eTransferQ").primitiveField();
        JJsigma.correctBoundaryConditions();
    }

    //- Form the current from the received E (see the block comment above)
    void eTransferProject
    (
        const Foam::fvMesh& mesh,
        Foam::electromagneticModel& em,
        const Foam::volVectorField& ERe,
        const Foam::volVectorField& EIm,
        const Foam::volScalarField& alpha1,
        const Foam::volVectorField& U,
        const Foam::volVectorField& Uexchange,
        const Foam::scalar floor,
        const Foam::word& mode,
        const bool newReference
    )
    {
        using namespace Foam;
        const dimensionSet dimE(1, 1, -3, 0, 0, -1, 0);   // V/m
        const dimensionedScalar sigma0(em.sigmaConst());
        const volScalarField a
        (
            "eTransferAlpha",
            min(max(alpha1, scalar(0)), scalar(1))
        );
        const volScalarField sigma("eTransferSigma", a*sigma0);
        const volScalarField sigmaS
        (
            "eTransferSigmaSolve",
            (floor + (1 - floor)*a)*sigma0
        );
        const volScalarField ratio("eTransferRatio", sigma/sigmaS);

        volScalarField& Q =
            eRegScalar(mesh, "eTransferQ", dimPower/dimVolume, false);
        Q.primitiveFieldRef() = 0;

        for (label part = 0; part < 2; part++)
        {
            const bool im = (part == 1);
            const word sfx(im ? "Im" : "Re");

            volVectorField& E =
                eRegVector(mesh, word("eTransferE") + sfx, dimE, true);
            E.primitiveFieldRef() = (im ? EIm : ERe).primitiveField();
            E.correctBoundaryConditions();

            const volVectorField& B = em.B(im);
            const volVectorField W("eTransferW", E + ((U - Uexchange) ^ B));

            surfaceScalarField psi
            (
                "eTransferPsi",
                fvc::interpolate(sigmaS*W) & mesh.Sf()
            );
            forAll(psi.boundaryField(), patchi)
            {
                if (!psi.boundaryField()[patchi].coupled())
                {
                    psi.boundaryFieldRef()[patchi] = 0;
                }
            }

            volScalarField& phi =
                eRegScalar(mesh, word("eTransferPhi") + sfx, dimE*dimLength, false);

            volVectorField Jnew("eTransferJ", sigma*W);
            scalar divRel = -1;

            if (mode == "none")
            {
                phi.primitiveFieldRef() = 0;
                phi.correctBoundaryConditions();
            }
            else
            {
                volScalarField rhs("eTransferRhs", fvc::div(psi));
                if (mode == "incremental")
                {
                    volScalarField& divRef =
                        eRegScalar
                        (
                            mesh, word("eTransferDivRef") + sfx, rhs.dimensions(), false
                        );
                    if (newReference)
                    {
                        divRef.primitiveFieldRef() = rhs.primitiveField();
                    }
                    rhs.primitiveFieldRef() -= divRef.primitiveField();
                }

                mesh.schemes().setFluxRequired(phi.name());
                fvScalarMatrix phiEqn
                (
                    fvm::laplacian(sigmaS, phi) == rhs
                );
                phiEqn.setReference(Pstream::master() ? 0 : -1, 0);
                phiEqn.solve();

                if (mode == "full")
                {
                    const surfaceScalarField Fl("eTransferF", psi - phiEqn.flux());
                    const surfaceVectorField Flv("eTransferFv", Fl*mesh.Cf());
                    Jnew =
                        ratio
                       *(
                            fvc::surfaceIntegrate(Flv)
                          - fvc::surfaceIntegrate(Fl)*mesh.C()
                        );
                    // discrete divergence left in the projected flux,
                    // relative to the largest face flux
                    const scalarField& Vc = mesh.V();
                    const scalarField netF
                    (
                        mag(fvc::surfaceIntegrate(Fl)().primitiveField())*Vc
                    );
                    divRel = gMax(netF)/max(gMax(mag(Fl.primitiveField())()), vSmall);
                }
                else // incremental
                {
                    const surfaceScalarField Fc("eTransferFc", -phiEqn.flux());
                    const surfaceVectorField Fcv("eTransferFcv", Fc*mesh.Cf());
                    Jnew +=
                        ratio
                       *(
                            fvc::surfaceIntegrate(Fcv)
                          - fvc::surfaceIntegrate(Fc)*mesh.C()
                        );
                }
            }

            // Joule heat sigma |E_total|^2 (no division by sigma)
            const volVectorField Etot("eTransferEtot", W - fvc::grad(phi));
            Q.primitiveFieldRef() +=
                0.5*sigma.primitiveField()*magSqr(Etot.primitiveField());

            // potential written in PotEre/PotEim for diagnostics
            volScalarField& PotE = em.PotE(im);
            PotE.primitiveFieldRef() = phi.primitiveField();

            volVectorField& J = em.J(im);
            J.primitiveFieldRef() = Jnew.primitiveField();
            J.correctBoundaryConditions();

            Info<< "eTransfer: " << mode << " " << sfx
                << ": max|phi| " << gMax(mag(phi.primitiveField())())
                << " V, max|E| melt "
                << gMax
                   (
                       pos0(a.primitiveField() - 0.999)
                      *mag(E.primitiveField())
                   )
                << " gas "
                << gMax
                   (
                       neg(a.primitiveField() - 0.01)
                      *mag(E.primitiveField())
                   )
                << " V/m, max|J| " << gMax(mag(J.primitiveField())())
                << ", net-flux residual " << divRel << endl;
        }
        Q.correctBoundaryConditions();
        eTransferForceAndHeat(mesh, em);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::incompressibleConductingVoFEtransfer::
incompressibleConductingVoFEtransfer(fvMesh& mesh)
:
    incompressibleConductingVoF(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::incompressibleConductingVoFEtransfer::
~incompressibleConductingVoFEtransfer()
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::word
Foam::solvers::incompressibleConductingVoFEtransfer::receivedField() const
{
    return electro.lookupOrDefault<word>("receivedField", "J");
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solvers::incompressibleConductingVoFEtransfer::solveElectromagnetics()
{
    // Electric-field transfer: the current is formed from the received E,
    // the current alpha1 and U (in fixCurrentsWithMixtureConductivity, cached
    // per time step); JxB and JJsigma are rewritten here because the
    // exchange-time electromagneticPredictor() recomputes them from J
    if (receivedField() == "E")
    {
        fixCurrentsWithMixtureConductivity();
        if (JReceivedRe_.valid() && mesh.foundObject<volScalarField>("eTransferQ"))
        {
            eTransferForceAndHeat(mesh, electro_);
        }
        return;
    }

    incompressibleConductingVoF::solveElectromagnetics();
}


void Foam::solvers::incompressibleConductingVoFEtransfer::
fixCurrentsWithMixtureConductivity()
{
    const word receivedField = this->receivedField();
    if (receivedField == "J")
    {
        incompressibleConductingVoF::fixCurrentsWithMixtureConductivity();
        return;
    }

    // As in incompressibleConductingVoF: before any exchange, J is left as is
    if (!JReceivedRe_.valid() || !alphaReceived_.valid())
    {
        return;
    }
    if (receivedField == "E")
    {
        if (!electro.isComplex() || !JReceivedIm_.valid())
        {
            FatalErrorInFunction
                << "receivedField E needs the harmonic model (E re and im)"
                << exit(FatalError);
        }
        const word mode =
            electro.lookupOrDefault<word>("eTransferProjection", "none");
        if (mode != "full" && mode != "none" && mode != "incremental")
        {
            FatalErrorInFunction
                << "eTransferProjection " << mode
                << " is not none, full or incremental" << exit(FatalError);
        }
        eTransferState& st = eState(mesh);
        const label timeIndex = runTime.timeIndex();
        if (timeIndex == st.doneTimeIndex && st.version == st.doneVersion)
        {
            return;
        }
        const bool newReference = (st.version != st.refVersion);
        eTransferProject
        (
            mesh, electro_, JReceivedRe_(), JReceivedIm_(), getAlpha1(), U_,
            U_old_, eTransferFloor, mode, newReference
        );
        st.refVersion = st.version;
        st.doneTimeIndex = timeIndex;
        st.doneVersion = st.version;
        return;
    }
    else
    {
        FatalErrorInFunction
            << "receivedField " << receivedField << " is not J or E"
            << exit(FatalError);
    }
}


void Foam::solvers::incompressibleConductingVoFEtransfer::storeReceivedCurrents
(
    bool imaginary
)
{
    if (receivedField() == "E")
    {
        // Electric-field transfer: a new received field invalidates the
        // cached current of this mesh
        eState(mesh).version++;
    }
    incompressibleConductingVoF::storeReceivedCurrents(imaginary);
}


// ************************************************************************* //
