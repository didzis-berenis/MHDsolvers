/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

#include "electromagneticModel.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvcSnGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
    defineTypeNameAndDebug(electromagneticModel, 0);
    defineRunTimeSelectionTable(electromagneticModel, fvMesh);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::volScalarField& Foam::electromagneticModel::lookupOrConstructScalar
(
    const fvMesh& mesh,
    const char* name,
    writeOption wo
)
{
    if (!mesh.objectRegistry::foundObject<volScalarField>(name))
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
                    IOobject::MUST_READ,
                    wo
                ),
                mesh
            )
        );

        // Transfer ownership of this object to the objectRegistry
        fPtr->store(fPtr);
    }

    return mesh.objectRegistry::lookupObjectRef<volScalarField>(name);
}

Foam::volVectorField& Foam::electromagneticModel::lookupOrConstructVector
(
    const fvMesh& mesh,
    const char* name,
    writeOption wo
)
{
    if (!mesh.objectRegistry::foundObject<volVectorField>(name))
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
                    IOobject::MUST_READ,
                    wo
                ),
                mesh
            )
        );

        // Transfer ownership of this object to the objectRegistry
        fPtr->store(fPtr);
    }

    return mesh.objectRegistry::lookupObjectRef<volVectorField>(name);
}

Foam::volScalarField& Foam::electromagneticModel::lookupOrConstructScalar
(
    const fvMesh& mesh,
    const char* name,
    dimensionedScalar value,
    readOption ro,
    writeOption wo
)
{
    if (!mesh.objectRegistry::foundObject<volScalarField>(name))
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
                    ro,
                    wo
                ),
                mesh,
                value
            )
        );

        // Transfer ownership of this object to the objectRegistry
        fPtr->store(fPtr);
    }

    return mesh.objectRegistry::lookupObjectRef<volScalarField>(name);
}

Foam::volVectorField& Foam::electromagneticModel::lookupOrConstructVector
(
    const fvMesh& mesh,
    const char* name,
    dimensionedVector value,
    readOption ro,
    writeOption wo
)
{
    if (!mesh.objectRegistry::foundObject<volVectorField>(name))
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
                    ro,
                    wo
                ),
                mesh,
                value
            )
        );

        // Transfer ownership of this object to the objectRegistry
        fPtr->store(fPtr);
    }

    return mesh.objectRegistry::lookupObjectRef<volVectorField>(name);
}

Foam::volVectorField& Foam::electromagneticModel::lookupOrConstructVector
(
    const fvMesh& mesh,
    const char* name,
    const volVectorField& field,
    readOption ro,
    writeOption wo
)
{
    if (!mesh.objectRegistry::foundObject<volVectorField>(name))
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
                    ro,
                    wo
                ),
                field
            )
        );

        // Transfer ownership of this object to the objectRegistry
        fPtr->store(fPtr);
    }

    return mesh.objectRegistry::lookupObjectRef<volVectorField>(name);
}

void Foam::electromagneticModel::constructScalar
(
    const fvMesh& mesh,
    const char* name
)
{
    if (!mesh.objectRegistry::foundObject<volScalarField>(name))
    {
        volScalarField* fPtr
        (
            new volScalarField
            (
                IOobject
                (
                    name,
                    mesh.time().name(),
                    mesh
                ),
                mesh
            )
        );

        // Transfer ownership of this object to the objectRegistry
        fPtr->store(fPtr);
    }
}

void Foam::electromagneticModel::constructVector
(
    const fvMesh& mesh,
    const char* name
)
{
    if (!mesh.objectRegistry::foundObject<volVectorField>(name))
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
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );

        // Transfer ownership of this object to the objectRegistry
        fPtr->store(fPtr);
    }
}


const Foam::electromagneticModel& Foam::electromagneticModel::lookupElectromagnetic
(
    const fvPatchScalarField& pf
)
{
    return pf.db().lookupObject<electromagneticModel>(typeName);
}

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

Foam::typeIOobject<Foam::IOdictionary>
Foam::electromagneticModel::readModelDict
(
    const objectRegistry& obr,
    const word& group,
    bool registerObject
)
{
    typeIOobject<IOdictionary> electromagneticDict
    (
        IOobject::groupName(typeName, group),
        obr.time().constant(),
        obr,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE,
        registerObject
    );
    return electromagneticDict;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::electromagneticModel::electromagneticModel
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    IOdictionary(readModelDict(mesh.thisDb(), phaseName, true)),
    mesh_(mesh),
    phaseName_(phaseName),
    JxB_(lookupOrConstructVector(mesh, "JxB")),
    JxB(JxB_),
    JJsigma_(lookupOrConstructScalar(mesh, "JJsigma")),
    JJsigma(JJsigma_),
    sigma_
    (
        lookupOrConstructScalar
        (
            mesh,
            "sigma",
            sigmaConst_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        )
    ),
    sigmaConst_
    (
        IOdictionary(readModelDict(mesh.thisDb(),phaseName)).found("sigma") ?
        dimensionedScalar
        (
            "sigma",
            pow3(dimTime)*dimCurrent*dimCurrent/dimMass/pow3(dimLength),
            IOdictionary(readModelDict(mesh.thisDb(),phaseName))
        ) :
        dimensionedScalar
        (
            "sigma",
            pow3(dimTime)*dimCurrent*dimCurrent/dimMass/pow3(dimLength),
            0
        )
    ),
    deltaU_
    (
        lookupOrConstructVector
        (
            mesh,
            "deltaU",
            dimensionedVector(dimVelocity,Foam::vector(0,0,0))
        )
    ),
    sigmaInv_
    (
        lookupOrConstructScalar
        (
            mesh,
            "sigmaInv",
            dimensionedScalar(dimMass*pow3(dimLength)/pow3(dimTime)/dimCurrent/dimCurrent,0)
        )
    )
{
    // Update sigma boundary fields
    if (sigmaConst_.value() > SMALL)
    {
        patchSigmaBoundaries();
    }
    // Update inverse of sigma
    forAll(sigma_, cellI)
    {
        sigmaInv_[cellI] = 
        sigma_[cellI] == 0 ? 
        0 : 1/sigma_[cellI];
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::electromagneticModel::updateDeltaU(volVectorField& Udiff)
{
    deltaU_ = Udiff;
}

void Foam::electromagneticModel::findDeltaJ(bool imaginary)
{
    /*---------------------------------------------------------------------------
    Correction to update electrical currents is based on the epotFoam solver,
    found in https://doi.org/10.13140/RG.2.2.12839.55201 (Chapter 4).
    ---------------------------------------------------------------------------*/
    //Interpolating cross product u x B over mesh faces
    surfaceScalarField psiUB = fvc::interpolate(deltaU_ ^ B(imaginary)) & mesh_.Sf();
    //Get reference for modification
    volScalarField& PotE = this->PotE(imaginary);
    //Sigma field
    volScalarField sigma_field(sigma_);
    //Poisson equation for electric potential
    fvScalarMatrix PotEEqn
    (
        fvm::laplacian(sigma_field,PotE)
        ==
        sigma_field*fvc::div(psiUB)
    );
    //Solving Poisson equation
    PotEEqn.solve();

    //Computation of current density at cell faces
    surfaceScalarField En = -(fvc::snGrad(PotE) * mesh_.magSf()) + psiUB;
    //Current density at face center
    surfaceVectorField Env = En * mesh_.Cf();

    //Get deltaJ reference (boundary conditions should come from J)
    volVectorField& JUB = this->deltaJ(imaginary);
    //Interpolation of current density at cell center
    JUB = sigma_field*(fvc::surfaceIntegrate(Env) - (fvc::surfaceIntegrate(En) * mesh_.C()) );
    //Update current density distribution and boundary conditions
    JUB.correctBoundaryConditions();
}

void Foam::electromagneticModel::predict()
{
    bool imaginary = isComplex();
    //Lorentz force term
    JxB_ =
    0.5*(
        (J() ^ B() )
        +(J(imaginary) ^ B(imaginary) )
    );
    JxB_.correctBoundaryConditions();
    //Joule heating
    //multiply by inverse of sigma to avoid division by zero
    JJsigma_ =
    0.5*(
        (J() & J())
        +(J(imaginary) & J(imaginary))
    )*sigmaInv();
    JJsigma_.correctBoundaryConditions();
}

void Foam::electromagneticModel::correct()
{
    bool imaginary = isComplex();
    //Get J difference by incorporating deltaU x B term
    findDeltaJ();
    if (imaginary)
    {
        findDeltaJ(imaginary);
    }
    volVectorField deltaJre = deltaJ();
    volVectorField deltaJim = deltaJ(imaginary);

    //Get references for modification
    //Lorentz force term
    JxB_ =
    0.5*(
        ((J()+deltaJre) ^ B() )
        +((J(imaginary)+deltaJim) ^ B(imaginary) )
    );
    JxB_.correctBoundaryConditions();
    //Joule heating
    //multiply by inverse of sigma to avoid division by zero
    JJsigma_ =
    0.5*(
        ((J()+deltaJre) & (J()+deltaJre))
        +((J(imaginary)+deltaJim) & (J(imaginary)+deltaJim))
    )*sigmaInv();
    JJsigma_.correctBoundaryConditions();
    //mark as corrected
    setCorrected();
}

bool Foam::electromagneticModel::read()
{
    return regIOobject::read();
}

Foam::volVectorField& Foam::electromagneticModel::getVectorFieldRef(const char* name)
{
    return mesh_.objectRegistry::lookupObjectRef<volVectorField>(name);
}

Foam::volScalarField& Foam::electromagneticModel::getScalarFieldRef(const char* name)
{
    return mesh_.objectRegistry::lookupObjectRef<volScalarField>(name);
}

void Foam::electromagneticModel::setCorrectElectromagnetics()
{
    correctElectromagnetics_ = true;
}

void Foam::electromagneticModel::setCorrected()
{
    correctElectromagnetics_ = false;
}

Foam::volScalarField& Foam::electromagneticModel::sigma()
{
    return sigma_;
}

void Foam::electromagneticModel::patchSigmaBoundaries()
{
    // Update sigma boundary fields.
    // Assign closest internal field value to boundary patch.
    // For interface walls (patches between regions),
    // sigma can have a different value from each side of the wall.
    volScalarField::Boundary& sigmaBf = sigma_.boundaryFieldRef();
    forAll(sigmaBf,patchi)
    {
        fvPatchScalarField& psigma = sigmaBf[patchi];
        psigma = psigma.patchInternalField();
    }
}

Foam::volScalarField& Foam::electromagneticModel::sigmaInv()
{
    return sigmaInv_;
}

// Read only access functions

const Foam::volVectorField& Foam::electromagneticModel::getVectorField(const char* name) const
{
    return mesh_.objectRegistry::lookupObject<volVectorField>(name);
}

const Foam::volScalarField& Foam::electromagneticModel::getScalarField(const char* name) const
{
    return mesh_.objectRegistry::lookupObject<volScalarField>(name);
}

bool Foam::electromagneticModel::correctElectromagnetics() const
{
    return correctElectromagnetics_;
}

Foam::tmp<Foam::scalarField> Foam::electromagneticModel::sigma(const label patchi) const
{
    return sigma_.boundaryField()[patchi].patchInternalField();
}

const Foam::volScalarField& Foam::electromagneticModel::sigmaInv() const
{
    return sigmaInv_;
}

const Foam::volScalarField& Foam::electromagneticModel::sigma() const
{
    return sigma_;
}

const Foam::dimensionedScalar Foam::electromagneticModel::sigmaConst() const
{
    return sigmaConst_;
}

// ************************************************************************* //
