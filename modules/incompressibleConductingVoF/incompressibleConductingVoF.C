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

#include "incompressibleConductingVoF.H"
#include "localEulerDdtScheme.H"
#include "fvCorrectPhi.H"
#include "geometricZeroField.H"
#include "findRefCell.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"
#include "fvcDiv.H"
#include "fvcSnGrad.H"
#include "surfaceInterpolate.H"
#include "surfaceTensionModel.H"
#include "contactAngleFvPatchScalarField.H"
#include "wedgeFvPatch.H"
#include "emptyFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(incompressibleConductingVoF, 0);
    addToRunTimeSelectionTable(solver, incompressibleConductingVoF, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::incompressibleConductingVoF::incompressibleConductingVoF(fvMesh& mesh)
:
    twoPhaseVoFSolver
    (
        mesh,
        autoPtr<twoPhaseVoFMixture>(new incompressibleTwoPhaseVoFMixture(mesh))
    ),

    mixture
    (
        refCast<incompressibleTwoPhaseVoFMixture>(twoPhaseVoFSolver::mixture)
    ),
    electroBase(mesh),

    p
    (
        IOobject
        (
            "p",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rho*buoyancy.gh
    ),

    alpha1_old_
    (
        IOobject
        (
            "alpha1_old",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha1
    ),

    alpha1_old(alpha1_old_),

    U_old_
    (
        IOobject
        (
            "U_old",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U_
    ),

    U_old(U_old_),

    curvatureSmoothing_
    (
        mixture.lookupOrDefault<label>("curvatureSmoothing", 0)
    ),

    interfaceK_
    (
        IOobject
        (
            "K",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless/dimLength, 0)
    ),

    pressureReference_
    (
        p,
        p_rgh,
        pimple.dict()
    ),

    momentumTransport
    (
        U,
        phi,
        alphaPhi1,
        mixture
    )
{
    label PotERefCell = 0;
    scalar PotERefValue = 0.0;
    setRefCell
    ( 
        electro.PotE(),
        pimple.dict(),
        PotERefCell,
        PotERefValue
    );

    if (electro.isComplex())
    {
        setRefCell
        ( 
            electro.PotE(true),
            pimple.dict(),
            PotERefCell,
            PotERefValue
        );
    }
    // Read the controls
    readControls();

    if (curvatureSmoothing_ < 0)
    {
        FatalErrorInFunction
            << "curvatureSmoothing " << curvatureSmoothing_
            << " in " << mixture.name() << " must be >= 0"
            << exit(FatalError);
    }
    if (curvatureSmoothing_ > 0)
    {
        // The smoothed curvature applies interfaceProperties'
        // contact-angle correction to its own face normal (updateInterfaceK)
        DynamicList<word> contactAnglePatches;
        forAll(alpha1.boundaryField(), patchi)
        {
            if
            (
                isA<contactAngleFvPatchScalarField>
                (
                    alpha1.boundaryField()[patchi]
                )
            )
            {
                contactAnglePatches.append(mesh.boundary()[patchi].name());
            }
        }
        if (contactAnglePatches.size())
        {
            Info<< "incompressibleConductingVoF: smoothed curvature corrects "
                << "the contact angle on the patches "
                << contactAnglePatches << " of " << alpha1.name() << endl;
        }
    }
    Info<< "incompressibleConductingVoF: curvatureSmoothing "
        << curvatureSmoothing_ << endl;
    updateInterfaceK();

    if (correctPhi || mesh.topoChanging())
    {
        rAU = new volScalarField
        (
            IOobject
            (
                "rAU",
                runTime.name(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimTime/dimDensity, 1)
        );
    }

    if (!runTime.restart() || !divergent())
    {
        correctUphiBCs(U_, phi_, true);

        fv::correctPhi
        (
            phi_,
            U,
            p_rgh,
            rAU,
            autoPtr<volScalarField>(),
            pressureReference(),
            pimple
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::incompressibleConductingVoF::~incompressibleConductingVoF()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::volScalarField& Foam::solvers::incompressibleConductingVoF::getAlpha1()
{
    return alpha1;
}

const Foam::volScalarField& Foam::solvers::incompressibleConductingVoF::getAlpha2()
{
    return alpha2;
}

const Foam::word& Foam::solvers::incompressibleConductingVoF::getPhase1Name()
{
    return mixture.phase1Name();
}

const Foam::word& Foam::solvers::incompressibleConductingVoF::getPhase2Name()
{
    return mixture.phase2Name();
}

void Foam::solvers::incompressibleConductingVoF::prePredictor()
{
    twoPhaseVoFSolver::prePredictor();

    const dimensionedScalar& rho1 = mixture.rho1();
    const dimensionedScalar& rho2 = mixture.rho2();

    // Calculate the mass-flux from the accumulated alphaPhi1
    rhoPhi = (alphaPhi1*(rho1 - rho2) + phi*rho2);

    if (pimple.predictTransport())
    {
        momentumTransport.predict();
    }
}


void Foam::solvers::incompressibleConductingVoF::pressureCorrector()
{
    incompressiblePressureCorrector(p);
}


void Foam::solvers::incompressibleConductingVoF::thermophysicalPredictor()
{}

const Foam::volScalarField& Foam::solvers::incompressibleConductingVoF::getAlpha1Old() const
{
    return alpha1_old_;
}

void Foam::solvers::incompressibleConductingVoF::storeAlpha1()
{
    alpha1_old_ = alpha1;
}

void Foam::solvers::incompressibleConductingVoF::storeU()
{
    U_old_ = U_;
}

const Foam::volScalarField& Foam::solvers::incompressibleConductingVoF::getNu()
{
    return mixture.nu();
}


void Foam::solvers::incompressibleConductingVoF::postCorrector()
{
    if (pimple.correctTransport())
    {
        momentumTransport.correct();
    }
}

void Foam::solvers::incompressibleConductingVoF::solveElectromagnetics()
{
    // Makes sure that currents stay in the conducting phase.
    fixCurrentsWithMixtureConductivity();

    Info<< nl << "Solving for electric potential" << endl;
    int nPote = 0;
    //Iterate until solved
    while (pimple.correct())
    {
        Info<< "Iteration " << ++nPote << endl;
        while (pimple.correctNonOrthogonal())
        {
            //Solve potential equation
            electro_.solve();
            if (pimple.finalNonOrthogonalIter() && electro.correctElectromagnetics())
            {
                // Update JxB and Joule heating terms with new potentials
                electro_.correct();
            }
        }
    }
    Info<< "Electric potential converged in " << nPote
    << " iterations" << endl;
}


void Foam::solvers::incompressibleConductingVoF::fixCurrentsWithMixtureConductivity()
{
    // The received current already carries the mixture conductivity
    // (alpha*sigma), so only rescale it for interface motion since the last
    // exchange. Rebuilding from the stored copy keeps repeated calls from
    // compounding the alpha weighting. Before any exchange, J is left as is.
    if (!JReceivedRe_.valid() || !alphaReceived_.valid())
    {
        return;
    }
    volScalarField alpha = getAlpha1();
    volScalarField a("emAlpha", min(max(alpha, scalar(0)), scalar(1)));
    // Capped at 1: a cell that was nearly gas at the exchange must not amplify J
    volScalarField ratio
    (
        "emAlphaRatio",
        min(a/max(alphaReceived_(), scalar(1e-3)), scalar(1))
    );

    volVectorField new_Jre = JReceivedRe_()*ratio;
    volVectorField& Jre = getJ();
    forAll(new_Jre, cellI)
    {
        Jre[cellI] = new_Jre[cellI];
    }
    Jre.correctBoundaryConditions();
    if (electro.isComplex() && JReceivedIm_.valid())
    {
        volVectorField new_Jim = JReceivedIm_()*ratio;
        volVectorField& Jim = getJ(true);
        forAll(new_Jim, cellI)
        {
            Jim[cellI] = new_Jim[cellI];
        }
        Jim.correctBoundaryConditions();
    }
}

void Foam::solvers::incompressibleConductingVoF::storeReceivedCurrents(bool imaginary)
{
    autoPtr<volVectorField>& JReceived = imaginary ? JReceivedIm_ : JReceivedRe_;
    JReceived.clear();
    JReceived.reset
    (
        new volVectorField
        (
            IOobject
            (
                imaginary ? "JReceivedIm" : "JReceivedRe",
                runTime.name(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            getJ(imaginary)
        )
    );
    if (!imaginary)
    {
        volScalarField alpha = getAlpha1();
        alphaReceived_.clear();
        alphaReceived_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "alphaReceived",
                    runTime.name(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                min(max(alpha, scalar(0)), scalar(1))
            )
        );
    }
}

Foam::tmp<Foam::volScalarField::Internal>
Foam::solvers::incompressibleConductingVoF::emForceLimiter()
{
    // Elmer's current does not follow alpha on this mesh near the interface,
    // and the small mixture density would turn it into a huge acceleration.
    // The limit comes from melt cells, where J does follow alpha. It uses the
    // unlimited force, so the limiter does not feed back on itself.
    const scalar meltAlpha = 0.9;
    const scalar gasAlphaCutoff = 1e-4;

    const volScalarField alpha(getAlpha1());
    const vectorField& JxB = electro.JxB();
    const scalarField magF(mag(JxB));
    const scalarField& rhoc = rho.primitiveField();

    scalar peak = 0;
    scalar stepPeak = 0;
    label hasReferenceCells = 0;
    forAll(magF, cellI)
    {
        if (alpha[cellI] >= gasAlphaCutoff)
        {
            stepPeak = max(stepPeak, magF[cellI]/rhoc[cellI]);
        }
        if (alpha[cellI] >= meltAlpha)
        {
            peak = max(peak, magF[cellI]/rhoc[cellI]);
            hasReferenceCells = 1;
        }
    }
    reduce(peak, maxOp<scalar>());
    reduce(stepPeak, maxOp<scalar>());
    reduce(hasReferenceCells, maxOp<label>());

    if (runTime.timeIndex() != emAccelerationTimeIndex_)
    {
        emAccelerationLimit_ = emAccelerationCurrent_;
        emAccelerationTimeIndex_ = runTime.timeIndex();
    }
    // Without melt reference cells, leave this step unbounded and use its
    // conducting-cell peak as the reference for the next step.
    emAccelerationCurrent_ = hasReferenceCells ? peak : stepPeak;
    const scalar limit =
        emAccelerationLimit_ > 0 ? emAccelerationLimit_ : peak;

    tmp<volScalarField::Internal> tweight
    (
        volScalarField::Internal::New
        (
            "emForceLimiter",
            mesh,
            dimensionedScalar(dimless, 1)
        )
    );
    scalarField& weight = tweight.ref();
    forAll(weight, cellI)
    {
        if (alpha[cellI] < gasAlphaCutoff)
        {
            weight[cellI] = 0;
        }
        else if
        (
            hasReferenceCells
         && magF[cellI] > rhoc[cellI]*limit
        )
        {
            weight[cellI] = rhoc[cellI]*limit/magF[cellI];
        }
    }
    return tweight;
}


void Foam::solvers::incompressibleConductingVoF::correctInterface()
{
    twoPhaseVoFSolver::correctInterface();
    updateInterfaceK();
}


void Foam::solvers::incompressibleConductingVoF::updateInterfaceK()
{
    if (curvatureSmoothing_ <= 0)
    {
        // The library curvature, unchanged (copy for output only)
        interfaceK_ =
            mesh.lookupObject<volScalarField>("interfaceProperties:K");
        return;
    }

    // Lafaurie et al. (1994) smoothing of a copy of alpha1: each pass
    // replaces the cell value by the area-weighted average of the linearly
    // interpolated face values. Wedge and empty patches are left out of the
    // average, so on an axisymmetric wedge the stencil is the meridional
    // one (the wedge faces would otherwise only add the cell's own value,
    // which near the axis would switch the smoothing off).
    volScalarField alphaS
    (
        IOobject
        (
            "alphaSmooth",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        alpha1
    );

    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const surfaceScalarField& magSf = mesh.magSf();

    for (label pass = 0; pass < curvatureSmoothing_; pass++)
    {
        const surfaceScalarField alphaSf(fvc::interpolate(alphaS));
        scalarField sumA(mesh.nCells(), Zero);
        scalarField sumW(mesh.nCells(), Zero);

        forAll(own, facei)
        {
            const scalar w = magSf[facei];
            const scalar wa = w*alphaSf[facei];
            sumA[own[facei]] += wa;
            sumW[own[facei]] += w;
            sumA[nei[facei]] += wa;
            sumW[nei[facei]] += w;
        }
        forAll(mesh.boundary(), patchi)
        {
            const fvPatch& fvp = mesh.boundary()[patchi];
            if (isA<wedgeFvPatch>(fvp) || isA<emptyFvPatch>(fvp))
            {
                continue;
            }
            const labelUList& fc = fvp.faceCells();
            const scalarField& pw = magSf.boundaryField()[patchi];
            const scalarField& pa = alphaSf.boundaryField()[patchi];
            forAll(fc, i)
            {
                sumA[fc[i]] += pw[i]*pa[i];
                sumW[fc[i]] += pw[i];
            }
        }
        alphaS.primitiveFieldRef() = sumA/sumW;
        alphaS.correctBoundaryConditions();
    }

    // interfaceProperties::calculateK with the smoothed alpha
    const volVectorField gradAlpha(fvc::grad(alphaS, "nHat"));
    const surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));
    surfaceVectorField nHatfv
    (
        gradAlphaf/(mag(gradAlphaf) + interface.deltaN())
    );

    // interfaceProperties::correctContactAngle on the smoothed face normal:
    // on every contactAngle patch of alpha1, reset nHatp to the contact angle
    // with the library's a12/b1/b2/det correction and the same deltaN, so
    // nHatf (below) and K carry the wall angle as at curvatureSmoothing 0.
    // The boundary condition of the real alpha1 (gradient, evaluate, alpha2)
    // is NOT set here: interface.correct() (twoPhaseVoFSolver::
    // correctInterface, called just before this function, and the
    // interfaceProperties constructor before the first call) has already
    // set it at this alpha1 from the unsmoothed normal. Setting it again
    // would evaluate the patch twice per call and change the alpha1 wall
    // values used by MULES and snGrad(alpha1).
    {
        const fvBoundaryMesh& boundary = mesh.boundary();
        surfaceVectorField::Boundary& nHatb = nHatfv.boundaryFieldRef();

        forAll(boundary, patchi)
        {
            if (isA<contactAngleFvPatchScalarField>(alpha1.boundaryField()[patchi]))
            {
                const contactAngleFvPatchScalarField& acap =
                    refCast<const contactAngleFvPatchScalarField>
                    (
                        alpha1.boundaryField()[patchi]
                    );

                fvsPatchVectorField& nHatp = nHatb[patchi];
                const scalarField cosTheta
                (
                    acap.cosTheta(U.boundaryField()[patchi], nHatp)
                );

                const vectorField nf(boundary[patchi].nf());

                const scalarField a12(nHatp & nf);
                const scalarField b1(cosTheta);

                scalarField b2(nHatp.size());
                forAll(b2, facei)
                {
                    b2[facei] = cos(acos(a12[facei]) - acos(cosTheta[facei]));
                }

                const scalarField det(1.0 - a12*a12);

                scalarField a((b1 - a12*b2)/det);
                scalarField b((b2 - a12*b1)/det);

                nHatp = a*nf + b*nHatp;
                nHatp /= (mag(nHatp) + interface.deltaN().value());
            }
        }
    }

    const surfaceScalarField nHatf(nHatfv & mesh.Sf());
    interfaceK_ = -fvc::div(nHatf);
}


Foam::tmp<Foam::surfaceScalarField>
Foam::solvers::incompressibleConductingVoF::surfaceTensionForce() const
{
    if (curvatureSmoothing_ <= 0)
    {
        return twoPhaseVoFSolver::surfaceTensionForce();
    }

    // As interfaceProperties::surfaceTensionForce, with the smoothed K and
    // the unsmoothed alpha1 in snGrad
    const surfaceTensionModel& sigmaModel =
        mesh.lookupObject<surfaceTensionModel>(surfaceTensionModel::typeName);

    return fvc::interpolate(sigmaModel.sigma()*interfaceK_)
       *fvc::snGrad(alpha1);
}

// ************************************************************************* //
