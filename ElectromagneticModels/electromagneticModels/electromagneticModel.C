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
//#include "fvModels.H"
//#include "fvConstraints.H"
//#include "findRefCell.H"
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
    regionRole_
    (
        this->lookupOrDefault<word>("regionRole", "")
    ),
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
        this->found("sigma") ?
        dimensionedScalar
        (
            "sigma",
            pow3(dimTime)*dimCurrent*dimCurrent/dimMass/pow3(dimLength),
            *this
        ) :
        dimensionedScalar
        (
            "sigma",
            pow3(dimTime)*dimCurrent*dimCurrent/dimMass/pow3(dimLength),
            0
        )
    ),
    /*deltaU_
    (
        lookupOrConstructVector
        (
            mesh,
            "deltaU",
            dimensionedVector(dimVelocity,Foam::vector(0,0,0))
        )
    ),*/
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
    // Update sigma boundary fields if necessary
    /*if (sigmaConst_.value() > SMALL)
    {
        patchSigmaBoundaries(sigma_);
    }*/
    // Update inverse of sigma
    forAll(sigma_, cellI)
    {
        sigmaInv_[cellI] = 
        sigma_[cellI] == 0 ? //Perhaps better to use "< SMALL" instead of "== 0"
        0 : 1/sigma_[cellI];
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::electromagneticModel::solve()
{
    findPotE();
    if (isComplex())
    {
        findPotE(true);
    }
    //Set as needs correction
    setCorrectElectromagnetics();
}
void Foam::electromagneticModel::findPotE(bool imaginary)
{
    /*---------------------------------------------------------------------------
    Correction to update electrical currents is based on the epotFoam solver,
    found in https://doi.org/10.13140/RG.2.2.12839.55201 (Chapter 4).
    ---------------------------------------------------------------------------*/
    //Interpolating cross product u x B over mesh faces
    surfaceScalarField psiUB = fvc::interpolate(deltaUxB(imaginary)) & mesh_.Sf();
    //Get reference for modification
    volScalarField& PotE = this->PotE(imaginary);
    //Poisson equation for electric potential
    fvScalarMatrix PotEEqn
    (
        fvm::laplacian(sigma_,PotE)
        ==
        sigma_*fvc::div(psiUB)
    );
    //Will add reference if needed
    label PotERefCell = 0;
    scalar PotERefValue = 0;
    PotEEqn.setReference(PotERefCell, PotERefValue);
    //Solving Poisson equation
    PotEEqn.solve();
}

void Foam::electromagneticModel::findDeltaJ(bool imaginary)
{
    // Add deltaU x B correction on top of externally calculated current density.
    // Is used in correct() method.
    /*---------------------------------------------------------------------------
    Correction to update electrical currents is based on the epotFoam solver,
    found in https://doi.org/10.13140/RG.2.2.12839.55201 (Chapter 4).
    ---------------------------------------------------------------------------*/
    //Get deltaJ reference (boundary conditions should come from J)
    volVectorField& JUB = this->deltaJ(imaginary);
    //Interpolating cross product u x B over mesh faces
    surfaceScalarField psiUB = fvc::interpolate(deltaUxB(imaginary)) & mesh_.Sf();//deltaU_ ^ B(imaginary)
    //Computation of current density at cell faces
    surfaceScalarField En = -(fvc::snGrad(PotE(imaginary)) * mesh_.magSf()) + psiUB;
    //Current density at face center
    surfaceVectorField Env = En * mesh_.Cf();

    //Interpolation of current density at cell center
    JUB = sigma_*(fvc::surfaceIntegrate(Env) - (fvc::surfaceIntegrate(En) * mesh_.C()) );
    //Update current density distribution and boundary conditions
    JUB.correctBoundaryConditions();
}

void Foam::electromagneticModel::findJ(bool imaginary)
{
    // Calculate current density internally.
    // To be used with predict() method.
    /*---------------------------------------------------------------------------
    Correction to update electrical currents is based on the epotFoam solver,
    found in https://doi.org/10.13140/RG.2.2.12839.55201 (Chapter 4).
    ---------------------------------------------------------------------------*/
    //Get J reference (boundary conditions should come from J)
    volVectorField& Jcorrect = this->J(imaginary);
    //Interpolating cross product u x B over mesh faces
    //surfaceScalarField psiUB = fvc::interpolate(UxB(imaginary)) & mesh_.Sf();//U_ ^ B(imaginary)
    //Computation of current density at cell faces
    surfaceScalarField En = -(fvc::snGrad(PotE(imaginary)) * mesh_.magSf());// + psiUB;
    //Current density at face center
    surfaceVectorField Env = En * mesh_.Cf();

    //Interpolation of current density at cell center
    Jcorrect = sigma_*(fvc::surfaceIntegrate(Env) - (fvc::surfaceIntegrate(En) * mesh_.C()) );
    //Update current density distribution and boundary conditions
    Jcorrect.correctBoundaryConditions();
}

void Foam::electromagneticModel::predict()
{
    // Updates Lorentz force and Joule heating
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
    // Calculates current correction and incorporates it in Lorentz force and Joule heating
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
    // Mark as corrected
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

/*void Foam::electromagneticModel::patchSigmaBoundaries(Foam::volScalarField& sigma)
{
    // Update sigma boundary fields.
    // Assign closest internal field value to boundary patch.
    // For interface walls (patches between regions),
    // sigma can have a different value from each side of the wall.
    volScalarField::Boundary& sigmaBf = sigma.boundaryFieldRef();
    forAll(sigmaBf,patchi)
    {
        fvPatchScalarField& psigma = sigmaBf[patchi];
        psigma = psigma.patchInternalField();
    }
}*/

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

Foam::word Foam::electromagneticModel::getRegionRole() const
{
    return regionRole_;
}

bool Foam::electromagneticModel::correctElectromagnetics() const
{
    return correctElectromagnetics_;
}

Foam::tmp<Foam::scalarField> Foam::electromagneticModel::sigma(const label patchi) const
{
    // Function patchInternalField() returns near-wall values mapped to the boundary.
    // In contrast, boundaryField()[patchi] will return the values at the boundary.
    return sigma_.boundaryField()[patchi].patchInternalField();
}

/*Foam::tmp<Foam::scalarField> Foam::electromagneticModel::sigmaInv(const label patchi) const
{
    //Return inverse of electric conductivity for patch
    return sigmaInv_.boundaryField()[patchi].patchInternalField();
}*/

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
