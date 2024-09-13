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
//#include "coupledCurrentDensityFvPatchVectorField.H"
#include "electromagneticModel.H"
#include "volFields.H"
#include "fvPatchFieldMapper.H"
#include "mappedPatchBase.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::coupledElectricPotentialFvPatchScalarField::getValues
(
    tmp<scalarField>& sigma,
    tmp<scalarField>& sigmaByDelta,
    tmp<scalarField>& EPot
) const
{
    const electromagneticModel& em =
        patch().boundaryMesh().mesh()
       .lookupType<electromagneticModel>();

    sigma = em.sigma(patch().index());
    sigmaByDelta = sigma()*patch().deltaCoeffs();
    EPot = patchInternalField();
}

void Foam::coupledElectricPotentialFvPatchScalarField::assign
(
    tmp<scalarField>& result,
    const tmp<scalarField>& field
) const
{
    if (result.valid())
    {
        result.clear();
    }

    if (field.isTmp())
    {
        result = field;
    }
    else
    {
        result = field().clone();
    }
}

Foam::word Foam::coupledElectricPotentialFvPatchScalarField::suffix() const
{
    const word ePotName = internalField().name();
    if ( ePotName != ePotName_ &&
        ePotName.size() >= ePotName_.size())
    {
        return ePotName.substr(ePotName_.size());
    }
    return "";
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
    mixedFvPatchScalarField(p, iF, dict, false)
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
        // default: valueFraction=1; (fixed value)
        valueFraction() = 1;
        refValue() = *this;
        refGrad() = 0;
    }
    //evaluate();// is called after all solvers have been constructed
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
    mixedFvPatchScalarField(psf, p, iF, mapper)
{}


Foam::coupledElectricPotentialFvPatchScalarField::
coupledElectricPotentialFvPatchScalarField
(
    const coupledElectricPotentialFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coupledElectricPotentialFvPatchScalarField::updateCoeffs()
{
    //See Foam::mixedFvPatchField<Type>::evaluate for details
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
        patchNbr.lookupPatchField<volScalarField, scalar>(ePotName_+suffix());

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

    // Get patch values
    tmp<scalarField> sigma;
    tmp<scalarField> sigmaByDelta;
    tmp<scalarField> EPot;
    getValues(sigma, sigmaByDelta, EPot);

    // Get patch values from deltaUxB field
    tmp<scalarField> deltaUxB;
    const fvPatchVectorField& deltaUxBp =
        patch().lookupPatchField<volVectorField, vector>(UxBname_+suffix());
    deltaUxB = deltaUxBp.patchInternalField() & patch().nf();

    // Get neighbour contributions
    tmp<scalarField> sigmaNbr;
    tmp<scalarField> sigmaByDeltaNbr;
    tmp<scalarField> EPotNbr;
    {
        tmp<scalarField> sigmaNbrPatch;
        tmp<scalarField> sigmaByDeltaNbrPatch;
        tmp<scalarField> EPotNbrPatch;
        coupledPotentialNbr.getValues(sigmaNbrPatch,sigmaByDeltaNbrPatch, EPotNbrPatch);
        assign(sigmaNbr, mpp.fromNeighbour(sigmaNbrPatch));
        assign(sigmaByDeltaNbr, mpp.fromNeighbour(sigmaByDeltaNbrPatch));
        assign(EPotNbr, mpp.fromNeighbour(EPotNbrPatch));
    }

    // Get neighbour contributions from deltaUxB field
    tmp<scalarField> deltaUxBNbr;
    const fvPatchVectorField& deltaUxBpNbr =
        patchNbr.lookupPatchField<volVectorField, vector>(UxBname_+suffix());
    {
        tmp<scalarField> deltaUxBNbrPatch;
        deltaUxBNbrPatch =
        deltaUxBpNbr.patchInternalField()
        &
        ( -deltaUxBpNbr.patch().nf() );//this->patch().nf() == -nbr.patch().nf()
        assign(deltaUxBNbr, mpp.fromNeighbour(deltaUxBNbrPatch));
    }

    this->refValue() =
    //sigma_nbr*ePot_nbr/delta_nbr + sigma*ePot/delta
    (sigmaByDeltaNbr()*EPotNbr()+sigmaByDelta()*EPot()+sigma()*deltaUxB()+sigmaNbr()*deltaUxBNbr()
    + SMALL*EPotNbr()//For the case of zero sigma in both regions
    )
    /
    //sigma_nbr/delta_nbr + sigma/delta
    (sigmaByDeltaNbr()+sigmaByDelta()+SMALL);

    mixedFvPatchScalarField::updateCoeffs();

    // Restore tag
    UPstream::msgType() = oldTag;
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
