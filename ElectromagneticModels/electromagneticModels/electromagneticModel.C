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
#include "zeroGradientFvPatchFields.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvcSnGrad.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(electromagneticModel, 0);
    //defineRunTimeSelectionTable(electromagneticModel, dictionary);

    Foam::autoPtr<electromagneticModel> electromagneticModel::New
    (
        const fvMesh& mesh,
        const word& phaseName
    )
    {
        const IOdictionary electromagneticDict
        (
            physicalProperties::findModelDict(mesh, phaseName)
        );
        const word modelType(electromagneticDict.lookup("electromagneticType"));

        Info<< "Selecting electromagnetics model " << modelType << endl;

        typename electromagneticModel::fvMeshConstructorTable::iterator
            cstrIter =
            electromagneticModel::fvMeshConstructorTablePtr_->find(modelType);

        if (cstrIter == electromagneticModel::fvMeshConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown " << electromagneticModel::typeName << " type "
                << modelType << nl << nl
                << "Valid " << electromagneticModel::typeName << " types are:" << nl
                << electromagneticModel::fvMeshConstructorTablePtr_->sortedToc() << nl
                << exit(FatalError);
        }

        return autoPtr<electromagneticModel>(cstrIter()(mesh, phaseName));
    }
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::volScalarField& Foam::electromagneticModel::lookupOrConstructScalar
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

    return mesh.objectRegistry::lookupObjectRef<volScalarField>(name);
}

Foam::volVectorField& Foam::electromagneticModel::lookupOrConstructVector
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
    return pf.db().lookupObject<electromagneticModel>(physicalProperties::typeName);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::electromagneticModel::electromagneticModel
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    mesh_(mesh),

    phaseName_(phaseName),

    physicalProperties(mesh, phaseName),

    JxB_(lookupOrConstructVector(mesh, "JxB")),

    JxB(JxB_),

    JJsigma_(lookupOrConstructScalar(mesh, "JJsigma")),

    JJsigma(JJsigma_),

    sigma_
    (
        IOobject
        (
            "sigma",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        sigmaConst_
    ),

    sigmaConst_
    (
        physicalProperties::lookupOrDefault<dimensionedScalar>
        (
            "sigma",
            dimensionedScalar
            (
                pow3(dimTime)*dimCurrent*dimCurrent/dimDensity,
                0
            )
        )
    ),

    sigmaInv_
    (
        IOobject
        (
            "sigmaInv",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimDensity/pow3(dimTime)/dimCurrent/dimCurrent,0)
    )
{
    //Update sigma patch fields
    if (sigmaConst_.value() > SMALL)
    {
        volScalarField::Boundary& sigmaBf =
            sigma_.boundaryFieldRef();
        forAll(sigmaBf,patchi)
        {
            fvPatchScalarField& psigma = sigmaBf[patchi];
            forAll(psigma, facei)
            {
                psigma[facei] = sigmaConst_.value();
            }
        }
    }

    forAll(sigma_, cellI)
    {
        sigmaInv_[cellI] = 
        sigma_[cellI] == 0 ? 
        0 : 1/sigma_[cellI];
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::electromagneticModel::updateDeltaU(const volVectorField& Udiff)
{
    deltaU.clear();
    deltaU = tmp<volVectorField>(Udiff);
}

Foam::volVectorField Foam::electromagneticModel::findDeltaJ(bool imaginary)
{
    //Interpolating cross product u x B over mesh faces
    surfaceScalarField psiUB = fvc::interpolate(deltaU ^ B(imaginary)) & mesh_.Sf();
    //Get reference for modification
    volScalarField& PotE = this->PotE(imaginary);
    //Const access
    volScalarField sigma(sigma_);
    //Poisson equation for electric potential
    fvScalarMatrix PotEEqn
    (
    fvm::laplacian(sigma,PotE)
    ==
    sigma*fvc::div(psiUB)
    );
    //Reference potential
    label PotERefCell = 0;
    scalar PotERefValue = 0.0;
    PotEEqn.setReference(PotERefCell, PotERefValue);
    //Solving Poisson equation
    PotEEqn.solve();

    //Computation of current density at cell faces
    surfaceScalarField En = -(fvc::snGrad(PotE) * mesh_.magSf()) + psiUB;
    //Current density at face center
    surfaceVectorField Env = En * mesh_.Cf();
    
    volVectorField JUB = deltaJ(imaginary).ref();
    //Interpolation of current density at cell center
    JUB = sigma*(fvc::surfaceIntegrate(Env) - (fvc::surfaceIntegrate(En) * mesh_.C()) );
    //Update current density distribution and boundary condition
    JUB.correctBoundaryConditions();

    return JUB;
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
    //Joule heating
    //multiply by inverse of sigma to avoid division by zero
    JJsigma_ =
    0.5*(
        (J() & J())
        +(J(imaginary) & J(imaginary))
    )*sigmaInv();
    setCorrected();
}

void Foam::electromagneticModel::correct()
{
    //volVectorField deltaU = U_-U_old_;
    //tmp<volVectorField> deltaU = deltaU()
    bool imaginary = isComplex();
    if (deltaU.valid())
    {
        volVectorField deltaJre = deltaJ();
        volVectorField deltaJim = deltaJ(imaginary);
        //Get J difference by incorporating deltaU x B term
        //volVectorField deltaJre = findDeltaJ(deltaU);
        deltaJre = findDeltaJ();
        //volVectorField deltaJim = deltaJre;
        deltaJim = deltaJre;
        if (imaginary)
        {
            deltaJim = findDeltaJ(imaginary);
        }

        //Get references for modification
        //Lorentz force term
        JxB_ =
        0.5*(
            ((J()+deltaJre) ^ B() )
            +((J(imaginary)+deltaJim) ^ B(imaginary) )
        );
        //Joule heating
        //multiply by inverse of sigma to avoid division by zero
        JJsigma_ =
        0.5*(
            ((J()+deltaJre) & (J()+deltaJre))
            +((J(imaginary)+deltaJim) & (J(imaginary)+deltaJim))
        )*sigmaInv();
        setCorrected();
    }
    else
    {
        predict();
    }
}

bool Foam::electromagneticModel::read()
{
    return regIOobject::read();
}

Foam::volVectorField& Foam::electromagneticModel::getVectorFromRegistry(const char* name)
{
    return mesh_.objectRegistry::lookupObjectRef<volVectorField>(name);
}

Foam::volScalarField& Foam::electromagneticModel::getScalarFromRegistry(const char* name)
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

bool Foam::electromagneticModel::correctElectromagnetics() const
{
    return correctElectromagnetics_;
}

Foam::volScalarField& Foam::electromagneticModel::sigma()
{
    return sigma_;
}

Foam::scalarField Foam::electromagneticModel::sigma(const label patchi) const
{
    return sigma_.boundaryField()[patchi];
}

Foam::volScalarField& Foam::electromagneticModel::sigmaInv()
{
    return sigmaInv_;
}
/*
Foam::volVectorField& Foam::electromagneticModel::JxB()
{
    return JxB_;
}

Foam::volScalarField& Foam::electromagneticModel::JJsigma()
{
    return JJsigma_;
}
*/
// Read only access functions

const Foam::tmp<Foam::volVectorField>& Foam::electromagneticModel::deltaJ(bool imaginary) const
{
    return deltaJ(imaginary);
}

const Foam::volScalarField& Foam::electromagneticModel::sigmaInv() const
{
    return sigmaInv();
}

const Foam::volScalarField& Foam::electromagneticModel::sigma() const
{
    return sigma();
}

const Foam::dimensionedScalar Foam::electromagneticModel::sigmaConst() const
{
    return sigmaConst_;
}

const Foam::volScalarField& Foam::electromagneticModel::PotE(bool imaginary) const
{
    return PotE(imaginary);
}

const Foam::volVectorField& Foam::electromagneticModel::J(bool imaginary) const
{
    return J(imaginary);
}

const Foam::volVectorField& Foam::electromagneticModel::B(bool imaginary) const
{
    return B(imaginary);
}
/*
const Foam::volVectorField& Foam::electromagneticModel::JxB() const
{
    return JxB();
}

const Foam::volScalarField& Foam::electromagneticModel::JJsigma() const
{
    return JJsigma();
}
*/

// ************************************************************************* //
