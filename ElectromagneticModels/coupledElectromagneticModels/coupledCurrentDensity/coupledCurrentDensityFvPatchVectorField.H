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
    Foam::coupledCurrentDensityFvPatchVectorField

Description

    Specifies gradient and temperature such that the equations are the same
    on both sides:
      - valueFraction = 1 -> fixed value
      - refValue = sigmaGradPotE,

    where sigmaGradPotE is the current density calculated from the electric
    potential gradient and the electric conductivity of the region.

SourceFiles
    coupledCurrentDensityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef coupledCurrentDensityFvPatchVectorField_H
#define coupledCurrentDensityFvPatchVectorField_H

#include "directionMixedFvPatchFields.H"
#include "coupledElectricPotentialFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class coupledCurrentDensityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class coupledCurrentDensityFvPatchVectorField
:
    public directionMixedFvPatchVectorField
{
    // Private Data

        //- Name of field on the neighbour region
        word ePotName_ = word::null;

        //- Coupled flag
        bool coupled_ = false;


protected:

    //- Get the patch sigma
    /*virtual void getValues
    (
        tmp<scalarField>& sigma
    ) const;*/

    //- Get current density from electric potential
    void getJfromPotential
    (
        const coupledElectricPotentialFvPatchScalarField& psf,
        tmp<scalarField>& sigmaGradPotE
    ) const;

public:

    //- Runtime type information
    TypeName("coupledCurrentDensity");


    // Constructors

        //- Construct as copy setting internal field reference
        coupledCurrentDensityFvPatchVectorField
        (
            const coupledCurrentDensityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );
        //- Construct from patch, internal field and dictionary
        coupledCurrentDensityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  coupledCurrentDensityFvPatchVectorField onto a
        //  new patch
        coupledCurrentDensityFvPatchVectorField
        (
            const coupledCurrentDensityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new coupledCurrentDensityFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        //- Initialize the coupled potential name
        virtual void initCoupledPotential();

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
