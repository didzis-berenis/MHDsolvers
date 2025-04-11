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

#include "electroBase.H"
#include "coupledCurrentDensityFvPatchVectorField.H"
#include "coupledElectricPotentialFvPatchScalarField.H"
#include "externalElectricPotentialFvPatchScalarField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(electroBase, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::electroBase::electroBase
(
    fvMesh& mesh
)
:
    electroPtr_(electromagneticModel::New(mesh)),
    electro_(electroPtr_()),
    electro(electroPtr_)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::electroBase::~electroBase()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::electroBase::solveElectromagnetics()
{
    //overridden in derived classes
}

void Foam::electroBase::electromagneticPredictor()
{
    electroPtr_->predict();
}

void Foam::electroBase::markCurrentReferenceAsSet()
{
    current_reference_set_ = true;
}

bool Foam::electroBase::currentReferenceSet()
{
    return current_reference_set_;
}

Foam::volVectorField& Foam::electroBase::getJ(bool imaginary)
{
    return electroPtr_->J(imaginary);
}

Foam::volVectorField& Foam::electroBase::getJref(bool imaginary)
{
    return electroPtr_->Jref(imaginary);
}

Foam::volVectorField& Foam::electroBase::getB(bool imaginary)
{
    return electroPtr_->B(imaginary);
}

void Foam::electroBase::initPotE(bool imaginary)
{
    volScalarField& PotE = electroPtr_->PotE(imaginary);
    volScalarField::Boundary& PotEBf = PotE.boundaryFieldRef();
    //Could also cast to mixedFvPatchScalarField,
    //which is the base class of coupledElectricPotentialFvPatchScalarField,
    //but this way it is more clear and excludes other uses of mixedFvPatchScalarField.
    forAll(PotEBf, patchi)
    {
        fvPatchScalarField& pPotE = PotEBf[patchi];
        if (isA<coupledElectricPotentialFvPatchScalarField>(pPotE) )
        {
            coupledElectricPotentialFvPatchScalarField& cpPotE =
            refCast<coupledElectricPotentialFvPatchScalarField>(pPotE);
            cpPotE.evaluate();
        }
    }
}


void Foam::electroBase::updatePotErefGrad(scalar newGrad, bool imaginary)
{
    volScalarField& PotE = electroPtr_->PotE(imaginary);
    volScalarField::Boundary& PotEBf = PotE.boundaryFieldRef();
    forAll(PotEBf, patchi)
    {
        fvPatchScalarField& pPotE = PotEBf[patchi];
        if (isA<coupledElectricPotentialFvPatchScalarField>(pPotE) )
        {
            coupledElectricPotentialFvPatchScalarField& cpPotE =
            refCast<coupledElectricPotentialFvPatchScalarField>(pPotE);
            if (cpPotE.getTerminalRole() == "terminal")
                cpPotE.refGrad() = newGrad;
        }
    }
}


void Foam::electroBase::updatePotErefValue(const word terminalName, scalar newValue, bool imaginary)
{
    volScalarField& PotE = electroPtr_->PotE(imaginary);
    volScalarField::Boundary& PotEBf = PotE.boundaryFieldRef();
    forAll(PotEBf, patchi)
    {
        fvPatchScalarField& pPotE = PotEBf[patchi];
        if (isA<coupledElectricPotentialFvPatchScalarField>(pPotE) )
        {
            coupledElectricPotentialFvPatchScalarField& cpPotE =
            refCast<coupledElectricPotentialFvPatchScalarField>(pPotE);
            if (cpPotE.getTerminalRole() == "terminal" && cpPotE.patch().name() == terminalName)
                cpPotE.refValue() = newValue;
        }
    }
}


Foam::scalar Foam::electroBase::getPotErefValue(const word terminalName, bool imaginary)
{
    volScalarField& PotE = electroPtr_->PotE(imaginary);
    volScalarField::Boundary& PotEBf = PotE.boundaryFieldRef();
    forAll(PotEBf, patchi)
    {
        fvPatchScalarField& pPotE = PotEBf[patchi];
        if (pPotE.patch().name() != terminalName) {continue;}
        if (isA<coupledElectricPotentialFvPatchScalarField>(pPotE) )
        {
            coupledElectricPotentialFvPatchScalarField& cpPotE =
            refCast<coupledElectricPotentialFvPatchScalarField>(pPotE);
            if (cpPotE.getTerminalRole() == "terminal")
                return gAverage(cpPotE.refValue());
        }
        if (isA<externalElectricPotentialFvPatchScalarField>(pPotE) )
        {
            externalElectricPotentialFvPatchScalarField& cpPotE =
            refCast<externalElectricPotentialFvPatchScalarField>(pPotE);
            if (cpPotE.hasPotE())
                return gAverage(cpPotE.refValue());
        }
    }
    if (hasBoundary(terminalName))
    {
        FatalIOError << "Boundary named " << terminalName << " is not defined as electric terminal!\n"
        << exit(FatalIOError);
    }
    else
    {
        FatalIOError << "Region doesn't have boundary named " << terminalName << "!\n"
        << exit(FatalIOError);
    }
    return -1;
}

bool Foam::electroBase::isElectricPotentialBoundary(const word terminalName)
{
    volScalarField& PotE = electroPtr_->PotE();
    volScalarField::Boundary& PotEBf = PotE.boundaryFieldRef();
    forAll(PotEBf, patchi)
    {
        fvPatchScalarField& pPotE = PotEBf[patchi];
        if (pPotE.patch().name() != terminalName)
            continue;
        if (isA<externalElectricPotentialFvPatchScalarField>(pPotE) )
        {
            externalElectricPotentialFvPatchScalarField& cpPotE =
            refCast<externalElectricPotentialFvPatchScalarField>(pPotE);
            if (cpPotE.hasPotE())
                return true;
        }
    }
    return false;
}

Foam::word Foam::electroBase::findGroundTerminal(const word terminalName)
{
    volScalarField& PotE = electroPtr_->PotE();
    volScalarField::Boundary& PotEBf = PotE.boundaryFieldRef();
    forAll(PotEBf, patchi)
    {
        fvPatchScalarField& pPotE = PotEBf[patchi];
        if (pPotE.patch().name() == terminalName)
            continue;
        if (isA<externalElectricPotentialFvPatchScalarField>(pPotE) )
        {
            externalElectricPotentialFvPatchScalarField& cpPotE =
            refCast<externalElectricPotentialFvPatchScalarField>(pPotE);
            if (cpPotE.hasPotE())
                return pPotE.patch().name();
        }
    }
    return "";
}

bool Foam::electroBase::hasBoundary(const word terminalName)
{
    volScalarField& PotE = electroPtr_->PotE();
    volScalarField::Boundary& PotEBf = PotE.boundaryFieldRef();
    forAll(PotEBf, patchi)
    {
        fvPatchScalarField& pPotE = PotEBf[patchi];
        if (pPotE.patch().name() == terminalName)
            return true;
        /*if (isA<coupledElectricPotentialFvPatchScalarField>(pPotE) )
        {
            coupledElectricPotentialFvPatchScalarField& cpPotE =
            refCast<coupledElectricPotentialFvPatchScalarField>(pPotE);
            if (cpPotE.getTerminalRole() == "terminal" && cpPotE.patch().name() == terminalName)
                return true;
        }*/
    }
    return false;
}


Foam::vector Foam::electroBase::getCenter(const word terminalName)
{
    volScalarField& PotE = electroPtr_->PotE();
    volScalarField::Boundary& PotEBf = PotE.boundaryFieldRef();
    forAll(PotEBf, patchi)
    {
        fvPatchScalarField& pPotE = PotEBf[patchi];
        if (pPotE.patch().name() == terminalName)
        {
            const scalarField surface = pPotE.patch().magSf();
            const vectorField points = pPotE.patch().Cf();
            const scalarField points_x = points & Foam::vector(1,0,0);
            const scalarField points_y = points & Foam::vector(0,1,0);
            const scalarField points_z = points & Foam::vector(0,0,1);
            scalar totVol = gSum(surface);
            scalar totVolx = gSum(surface*points_x)/totVol;
            scalar totVoly = gSum(surface*points_y)/totVol;
            scalar totVolz = gSum(surface*points_z)/totVol;
            Foam::vector centerPosition(totVolx,totVoly,totVolz);
            return centerPosition;
        }
    }
    FatalIOError << "Region doesn't have boundary named " << terminalName << "!\n"
    << exit(FatalIOError);
    return Foam::vector(0,0,0);
}


void Foam::electroBase::initDeltaJ(bool imaginary)
{
        volVectorField& deltaJ = electroPtr_->deltaJ(imaginary);
        volVectorField::Boundary& deltaJBf = deltaJ.boundaryFieldRef();
        forAll(deltaJBf, patchi)
        {
            fvPatchVectorField& pDeltaJ = deltaJBf[patchi];
            if (isA<coupledCurrentDensityFvPatchVectorField>(pDeltaJ) )
            {
                //derived from directionMixedFvPatchVectorField
                coupledCurrentDensityFvPatchVectorField& cpDeltaJ =
                refCast<coupledCurrentDensityFvPatchVectorField>(pDeltaJ);
                //Switch to coupling mode
                cpDeltaJ.initCoupling();
                cpDeltaJ.evaluate();
            }
        }
}

void Foam::electroBase::initJ(bool imaginary)
{
        volVectorField& J = electroPtr_->J(imaginary);
        volVectorField::Boundary& JBf = J.boundaryFieldRef();
        forAll(JBf, patchi)
        {
            fvPatchVectorField& pJ = JBf[patchi];
            if (isA<coupledCurrentDensityFvPatchVectorField>(pJ) )
            {
                //derived from directionMixedFvPatchVectorField
                coupledCurrentDensityFvPatchVectorField& cpJ =
                refCast<coupledCurrentDensityFvPatchVectorField>(pJ);
                //Switch to coupling mode
                cpJ.initCoupling();
                cpJ.evaluate();
            }
        }
}

// ************************************************************************* //
