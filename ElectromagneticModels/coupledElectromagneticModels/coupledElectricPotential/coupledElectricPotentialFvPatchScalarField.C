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

#include "coupledElectricPotentialFvPatchScalarField.H"
#include "electromagneticModel.H"
#include "volFields.H"
#include "fvPatchFieldMapper.H"
#include "mappedPatchBase.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::coupledElectricPotentialFvPatchScalarField::getThis
(
    tmp<scalarField>& sigma,
    tmp<scalarField>& ePotByDelta
) const
{
    const electromagneticModel& em =
        patch().boundaryMesh().mesh()
       .lookupType<electromagneticModel>();

    sigma = em.sigma(patch().index());

    ePotByDelta = patchInternalField()*patch().deltaCoeffs();
}


void Foam::coupledElectricPotentialFvPatchScalarField::getNbr
(
    tmp<scalarField>& sigmaByDeltaNbr,
    tmp<scalarField>& sigmaEPotByDeltaNbr
) const
{
    const electromagneticModel& em =
        patch().boundaryMesh().mesh()
       .lookupType<electromagneticModel>();

    sigmaByDeltaNbr = em.sigma(patch().index())*patch().deltaCoeffs();

    sigmaEPotByDeltaNbr = sigmaByDeltaNbr()*patchInternalField();
}


void Foam::coupledElectricPotentialFvPatchScalarField::add
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

Foam::coupledElectricPotentialFvPatchScalarField::
coupledElectricPotentialFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict, false),
    //JName_(dict.lookupOrDefault<word>("J", "J")),
    ePotnbrName_(dict.lookupOrDefault<word>("PotE", "PotE"))
{
    mappedPatchBase::validateMapForField
    (
        *this,
        iF,
        dict,
        mappedPatchBase::from::differentPatch
    );

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0;
        valueFraction() = 1;
    }
}


Foam::coupledElectricPotentialFvPatchScalarField::
coupledElectricPotentialFvPatchScalarField
(
    const coupledElectricPotentialFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper),
   //JName_(psf.JName_),
    ePotnbrName_(psf.ePotnbrName_)
{}


Foam::coupledElectricPotentialFvPatchScalarField::
coupledElectricPotentialFvPatchScalarField
(
    const coupledElectricPotentialFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF),
    //JName_(psf.JName_),
    ePotnbrName_(psf.ePotnbrName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coupledElectricPotentialFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
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
        patchNbr.lookupPatchField<volScalarField, scalar>(ePotnbrName_);

    if (!isA<coupledElectricPotentialFvPatchScalarField>(ePotpNbr))
    {
        FatalErrorInFunction
            << "Patch field for " << internalField().name() << " on "
            << this->patch().name() << " is of type "
            << coupledElectricPotentialFvPatchScalarField::typeName
            << endl << "The neighbouring patch field "
            << internalField().name() << " on "
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
    tmp<scalarField> ePotByDelta;
    // Get patch values
    getThis(sigma, ePotByDelta);

    tmp<scalarField> sigmaByDelta;
    tmp<scalarField> sigmaEPotByDelta;
    // Add neighbour contributions
    {
        tmp<scalarField> sigmaByDeltaNbr;
        tmp<scalarField> sigmaEPotByDeltaNbr;
        coupledPotentialNbr.getNbr(sigmaByDeltaNbr, sigmaEPotByDeltaNbr);
        //
        add(sigmaEPotByDelta, mpp.fromNeighbour(sigmaEPotByDeltaNbr));
        add(sigmaByDelta, mpp.fromNeighbour(sigmaByDeltaNbr));
    }

    //const Field<vector>& JpNbr =
    //    patchNbr.lookupPatchField<volVectorField, vector>(JName_);
    //const Field<vector> nfNbr(mpp.nbrPolyPatch().nf());
    //const scalarField nJpNbr(JpNbr & nfNbr);

    // default: 0; if sigma->0 => valueFraction=1
    this->valueFraction() = 1 - sigma()/(sigma() + SMALL);
    // if sigma->: ePot = ePotNbr; if sigmaNbr-> => ePot = 0
    this->refValue() = sigmaEPotByDelta()/(sigmaByDelta()+SMALL);
    // default: using gradient grad(ePot) = grad(ePotNbr)*sigmaNbr/sigma
    // if sigmaNbr-> => grad(ePot) = 0
    this->refGrad() = sigmaEPotByDelta()/(sigma()+SMALL);

    mixedFvPatchScalarField::updateCoeffs();

    // Restore tag
    UPstream::msgType() = oldTag;
}


void Foam::coupledElectricPotentialFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "PotE", "PotE", ePotnbrName_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        coupledElectricPotentialFvPatchScalarField
    );
}


// ************************************************************************* //
