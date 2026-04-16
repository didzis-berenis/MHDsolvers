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

#include "uxbElectricPotentialFvPatchScalarField.H"
#include "electromagneticModel.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

Foam::word Foam::uxbElectricPotentialFvPatchScalarField::suffix() const
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

Foam::uxbElectricPotentialFvPatchScalarField::
uxbElectricPotentialFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict, false),
    UxBname_(dict.lookupOrDefault<word>("UxBname","deltaUxB"))
{
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume zeroGradient.
        refValue() = 0;
        refGrad() = 0;
        valueFraction() = 0;
    }
}


Foam::uxbElectricPotentialFvPatchScalarField::
uxbElectricPotentialFvPatchScalarField
(
    const uxbElectricPotentialFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::uxbElectricPotentialFvPatchScalarField::
uxbElectricPotentialFvPatchScalarField
(
    const uxbElectricPotentialFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::uxbElectricPotentialFvPatchScalarField::map
(
    const fvPatchScalarField& ptf,
    const fvPatchFieldMapper& mapper
)
{
    mixedFvPatchScalarField::map(ptf, mapper);

    const uxbElectricPotentialFvPatchScalarField& tiptf =
        refCast<const uxbElectricPotentialFvPatchScalarField>(ptf);
}


void Foam::uxbElectricPotentialFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    mixedFvPatchScalarField::reset(ptf);

    const uxbElectricPotentialFvPatchScalarField& tiptf =
        refCast<const uxbElectricPotentialFvPatchScalarField>(ptf);
}

void Foam::uxbElectricPotentialFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalarField& ePotP(*this);
    scalarField gradPotE(ePotP.size(), 0);

    // Get UxB from closest value (since the boundary conditions are not set for this field)
    word UxBname = UxBname_ == word("deltaUxB") ? word(UxBname_+suffix()) : UxBname_;
    const fvPatchVectorField& deltaUxBp =
        patch().lookupPatchField<volVectorField, vector>(UxBname);
    tmp<scalarField> deltaUxB;
    deltaUxB = deltaUxBp.patchInternalField() & patch().nf();

    this->refGrad() = deltaUxB();
    
    mixedFvPatchScalarField::updateCoeffs();

}


void Foam::uxbElectricPotentialFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);

    writeEntry(os, "UxBname", UxBname_ == word("deltaUxB") ? word(UxBname_+suffix()) : UxBname_);
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
        uxbElectricPotentialFvPatchScalarField
    );
}


// ************************************************************************* //
