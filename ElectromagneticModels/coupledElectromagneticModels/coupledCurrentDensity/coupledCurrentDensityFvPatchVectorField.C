/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "coupledCurrentDensityFvPatchVectorField.H"
#include "coupledElectricPotentialFvPatchScalarField.H"
#include "electromagneticModel.H"
#include "volFields.H"
#include "fvPatchFieldMapper.H"
#include "mappedPatchBase.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::coupledCurrentDensityFvPatchVectorField::getThis
(
    tmp<scalarField>& sigma//,
    //tmp<scalarField>& ePotByDelta
) const
{
    const electromagneticModel& em =
        patch().boundaryMesh().mesh()
       .lookupType<electromagneticModel>();

    sigma = em.sigma(patch().index());

    //ePotByDelta = patchInternalField()*patch().deltaCoeffs();
}

const Foam::word Foam::coupledCurrentDensityFvPatchVectorField::getCoupledPotentialName() const
{
    const electromagneticModel& em =
        patch().boundaryMesh().mesh()
       .lookupType<electromagneticModel>();

    return em.getCoupledPotentialName(internalField().name());
}

void Foam::coupledCurrentDensityFvPatchVectorField::add
(
    tmp<scalarField>& result,
    const tmp<scalarField>& field
) const
{
    if (result.valid())
    {
        result.ref() += field;
    }
    else
    {
        if (field.isTmp())
        {
            result = field;
        }
        else
        {
            result = field().clone();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coupledCurrentDensityFvPatchVectorField::
coupledCurrentDensityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    directionMixedFvPatchVectorField(p, iF),
    //JName_(dict.lookupOrDefault<word>("J", "J"))
    ePotName_(getCoupledPotentialName())
{
    Info << "initializing coupledCurrentDensity for patch " << this->patch().name() << " and field " << internalField().name() << endl;
    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = vectorField("refValue", dict, p.size());
        refGrad() = vectorField("refGradient", dict, p.size());
        valueFraction() = symmTensorField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        //refValue() = *this;
        refValue() = Zero;//normalCurrent();
        //valueFraction() is symmTensorField:
        //"xx", "xy", "xz",
        //"yy", "yz",
        //      "zz"
        valueFraction() = sqr(patch().nf());//fraction == normal direction
        refGrad() = Zero;
    }
    //evaluate() is called after all solvers have been constructed
}


Foam::coupledCurrentDensityFvPatchVectorField::
coupledCurrentDensityFvPatchVectorField
(
    const coupledCurrentDensityFvPatchVectorField& psf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directionMixedFvPatchVectorField(psf, p, iF, mapper),
    //JName_(psf.JName_)
    ePotName_(psf.ePotName_)
{}


Foam::coupledCurrentDensityFvPatchVectorField::
coupledCurrentDensityFvPatchVectorField
(
    const coupledCurrentDensityFvPatchVectorField& psf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(psf, iF),
    //JName_(psf.JName_)
    ePotName_(psf.ePotName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coupledCurrentDensityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    Info << "updating coupledCurrentDensity for patch " << this->patch().name() << " and field " << internalField().name() << endl;
    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag + 1;

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp = mappedPatchBase::getMap(patch().patch());
    const label patchiNbr = mpp.nbrPolyPatch().index();
    const fvPatch& patchNbr =
        refCast<const fvMesh>(mpp.nbrMesh()).boundary()[patchiNbr];

    const fvPatchScalarField& ePotpNbr =
        patchNbr.lookupPatchField<volScalarField, scalar>(ePotName_);

    if (!isA<coupledElectricPotentialFvPatchScalarField>(ePotpNbr))
    {
        FatalErrorInFunction
            << "Patch field for " << ePotName_ << " on "
            << this->patch().name() << " is of type "
            << coupledElectricPotentialFvPatchScalarField::typeName
            << endl << "The neighbouring patch field "
            << ePotName_ << " on "
            << patchNbr.name() << " is required to be the same, but is "
            << "currently of type " << ePotpNbr.type() << exit(FatalError);
    }

    const coupledElectricPotentialFvPatchScalarField& coupledPotentialNbr =
        refCast<const coupledElectricPotentialFvPatchScalarField>(ePotpNbr);

    //if (Js_.valid())
    //{
    //    sumJ += Js_();
    //}
    //const Field<vector>& Jp =
    //    patch().template lookupPatchField<volVectorField, vector>(JName_);
    //const Field<vector> nf(patch().nf());
    //const scalarField nJp(Jp & nf);

    tmp<scalarField> sigma;
    //tmp<scalarField> ePotByDelta;
    // Get patch values
    getThis(sigma);//, ePotByDelta);

    //tmp<scalarField> sigmaByDelta;
    tmp<scalarField> sigmaEPotByDelta;
    // Add neighbour contributions
    {
        tmp<scalarField> sigmaEPotByDeltaNbr;
        sigmaEPotByDeltaNbr =
        coupledPotentialNbr.getNbr();//sigmaByDeltaNbr, 
        //sigmaEPotByDeltaNbr);
        //
        add(sigmaEPotByDelta, mpp.fromNeighbour(sigmaEPotByDeltaNbr));
        //add(sigmaByDelta, mpp.fromNeighbour(sigmaByDeltaNbr));
    }

    //const Field<vector>& JpNbr =
    //    patchNbr.lookupPatchField<volVectorField, vector>(JName_);
    //const Field<vector> nfNbr(mpp.nbrPolyPatch().nf());
    //const scalarField nJpNbr(JpNbr & nfNbr);

    // default: J . n = Jnbr . n = grad(ePotNbr)*sigmaNbr . n
    // if sigmaNbr->0 => Jnbr . n = 0
    // if sigma->0 => sigma()/(sigma() + SMALL) = 0 => J . n = 0
    this->refValue() = (sigma()/(sigma() + SMALL)) * sigmaEPotByDelta() * patch().nf();

    directionMixedFvPatchVectorField::updateCoeffs();

    // Restore tag
    UPstream::msgType() = oldTag;
}


void Foam::coupledCurrentDensityFvPatchVectorField::write
(
    Ostream& os
) const
{
    directionMixedFvPatchVectorField::write(os);
    writeEntryIfDifferent<word>(os, "PotE", "PotE", ePotName_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        coupledCurrentDensityFvPatchVectorField
    );
}


// ************************************************************************* //
