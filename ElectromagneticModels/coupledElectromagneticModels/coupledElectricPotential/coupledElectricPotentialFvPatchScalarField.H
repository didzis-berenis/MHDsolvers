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

Class
    Foam::coupledElectricPotentialFvPatchScalarField

Description
    Boundary condition for electrical potential, to be used for electric contact
    with another region in an electromagnetic case.

    Specifies potential such that the equations are the same
    on both sides:
      - refGradient = 0
      - refValue = neighbour value
      - valueFraction = 1 (fixed value)
Usage

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            coupledPotential;
    }
    \endverbatim

    Needs to be on underlying mapped(Wall)FvPatch.

    The patch electric conductivity sigma is obtained from the region
    electromagneticModel so that this boundary condition can be applied
    directly to either fluid or solid regions.

SourceFiles
    coupledElectricPotentialFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef coupledElectricPotentialFvPatchScalarField_H
#define coupledElectricPotentialFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class coupledElectricPotentialFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class coupledElectricPotentialFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private Data

        //- Name of field on the neighbour region
        const word ePotName_;


protected:

    //- Get the neighbour patch sigma and sigma*PotE
    virtual void getNbr
    (
        tmp<scalarField>& sigmaNbr,
        tmp<scalarField>& sigmaEPotNbr
    ) const;

    //- Add field to result which may have not been previously set
    void add
    (
        tmp<scalarField>& result,
        const tmp<scalarField>& field
    ) const;


public:

    //- Runtime type information
    TypeName("coupledPotential");

    // Constructors

        //- Construct as copy setting internal field reference
        coupledElectricPotentialFvPatchScalarField
        (
            const coupledElectricPotentialFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );
        //- Construct from patch, internal field and dictionary
        coupledElectricPotentialFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  coupledElectricPotentialFvPatchScalarField onto a
        //  new patch
        coupledElectricPotentialFvPatchScalarField
        (
            const coupledElectricPotentialFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new coupledElectricPotentialFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Return the field from the object registry
        virtual const scalarField& getNbr() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
