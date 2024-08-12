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

#include "externalElectricPotentialFvPatchScalarField.H"
#include "electromagneticModel.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::externalElectricPotentialFvPatchScalarField::
externalElectricPotentialFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict, false),
    haveI_(dict.found("I")),
    I_(haveI_ ? dict.lookup<scalar>("I") : NaN),
    haveJ_(dict.found("J")),
    J_(haveJ_ ? scalarField("J", dict, p.size()) : scalarField()),
    ePotExt_(havePot_ ? Function1<scalar>::New("PotEa", dict).ptr() : nullptr)
    /*,
    thicknessLayers_
    (
        dict.lookupOrDefault<scalarList>("thicknessLayers", scalarList())
    ),
    kappaLayers_
    (
        dict.lookupOrDefault<scalarList>("kappaLayers", scalarList())
    )*/
{
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (!haveI_ && !haveJ_ && !havePot_)
    {
        FatalIOErrorInFunction(dict)
            << "One or more of Q (heat power), q (heat flux), and h (heat "
            << "transfer coefficient) must be specified"
            << exit(FatalIOError);
    }
/*
    if (thicknessLayers_.size() != kappaLayers_.size())
    {
        FatalIOErrorInFunction(dict)
            << "If either thicknessLayers or kappaLayers is specified, then "
            << "both must be specified and be lists of the same length "
            << exit(FatalIOError);
    }
*/
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


Foam::externalElectricPotentialFvPatchScalarField::
externalElectricPotentialFvPatchScalarField
(
    const externalElectricPotentialFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    haveI_(ptf.haveI_),
    I_(ptf.I_),
    haveJ_(ptf.haveJ_),
    J_(haveJ_ ? mapper(ptf.J_)() : scalarField()),
    havePot_(ptf.havePot_),
    ePotExt_(ptf.ePotExt_, false)
    //thicknessLayers_(ptf.thicknessLayers_),
    //kappaLayers_(ptf.kappaLayers_)
{}


Foam::externalElectricPotentialFvPatchScalarField::
externalElectricPotentialFvPatchScalarField
(
    const externalElectricPotentialFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF),
    haveI_(tppsf.haveI_),
    I_(tppsf.I_),
    haveJ_(tppsf.haveJ_),
    J_(tppsf.J_),
    havePot_(tppsf.havePot_),
    ePotExt_(tppsf.ePotExt_, false)
    //thicknessLayers_(tppsf.thicknessLayers_),
    //kappaLayers_(tppsf.kappaLayers_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::externalElectricPotentialFvPatchScalarField::map
(
    const fvPatchScalarField& ptf,
    const fvPatchFieldMapper& mapper
)
{
    mixedFvPatchScalarField::map(ptf, mapper);

    const externalElectricPotentialFvPatchScalarField& tiptf =
        refCast<const externalElectricPotentialFvPatchScalarField>(ptf);

    if (haveJ_)
    {
        mapper(J_, tiptf.J_);
    }
/*
    if (havePot_)
    {
        mapper(ePotExt_, tiptf.ePotExt_);
    }
*/
}


void Foam::externalElectricPotentialFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    mixedFvPatchScalarField::reset(ptf);

    const externalElectricPotentialFvPatchScalarField& tiptf =
        refCast<const externalElectricPotentialFvPatchScalarField>(ptf);

    if (haveJ_)
    {
        J_.reset(tiptf.J_);
    }
/*
    if (havePot_)
    {
        ePotExt_.reset(tiptf.ePotExt_);
    }
*/
}


void Foam::externalElectricPotentialFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalarField& ePotP(*this);

    const electromagneticModel& em =
        patch().boundaryMesh().mesh()
       .lookupType<electromagneticModel>();

    const scalarField sigma(em.sigma(patch().index()));
    const scalar ePotExt = ePotExt_->value(this->db().time().userTimeValue());

    // Compute the total non-convective heat flux
    scalarField ePotTot(ePotP.size(), 0);
    if (haveI_)
    {
        ePotTot += sigma*I_/gSum(patch().magSf());
    }
    if (haveJ_)
    {
        ePotTot += sigma*J_;
    }
    if (!havePot_)
    {
        ePotTot += ePotExt;
    }

    // Evaluate
    //if (!haveh_)
    //{
        refGrad() = ePotTot/sigma;
        refValue() = ePotP;
        valueFraction() = 0;
    //}
/*
    else
    {
        scalar totalSolidRes = 0;
        if (thicknessLayers_.size())
        {
            forAll(thicknessLayers_, iLayer)
            {
                const scalar l = thicknessLayers_[iLayer];
                if (kappaLayers_[iLayer] > 0)
                {
                    totalSolidRes += l/kappaLayers_[iLayer];
                }
            }
        }

        const scalar Ta = Ta_->value(this->db().time().userTimeValue());

        const scalarField hp
        (
            h_ + totalSolidRes
        );

        const scalarField hpTa(hp*Ta);

        const scalarField kappaDeltaCoeffs
        (
            kappa*patch().deltaCoeffs()
        );

        refGrad() = 0;
        forAll(Tp, i)
        {
            if (qTot[i] < 0)
            {
                const scalar hpmqTot = hp[i] - qTot[i]/Tp[i];
                refValue()[i] = hpTa[i]/hpmqTot;
                valueFraction()[i] = hpmqTot/(hpmqTot + kappaDeltaCoeffs[i]);
            }
            else
            {
                refValue()[i] = (hpTa[i] + qTot[i])/hp[i];
                valueFraction()[i] = hp[i]/(hp[i] + kappaDeltaCoeffs[i]);
            }
        }
    }
*/

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        const scalar Q = gSum(sigma*patch().magSf()*snGrad());

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " :"
            << " heat transfer rate:" << Q
            << " walltemperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }
}


void Foam::externalElectricPotentialFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);

    if (haveI_)
    {
        writeEntry(os, "I", I_);
    }

    if (haveJ_)
    {
        writeEntry(os, "J", J_);
    }

    if (havePot_)
    {
        writeEntry(os, ePotExt_());
/*
        writeEntry(os, "h", h_);
        writeEntryIfDifferent
        (
            os,
            "thicknessLayers",
            scalarList(),
            thicknessLayers_
        );
        writeEntryIfDifferent
        (
            os,
            "kappaLayers",
            scalarList(),
            kappaLayers_
        );
*/
    }
    writeEntry(os, "refValue", refValue());
    writeEntry(os, "refGradient", refGrad());
    writeEntry(os, "valueFraction", valueFraction());
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        externalElectricPotentialFvPatchScalarField
    );
}


// ************************************************************************* //
