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
      - refValue = calculated as patch_value from
        sigma_1*( -grad(PotE_1) + UxB_1 ) = sigma_2*( -grad(PotE_2) + UxB_2 ) where sigma_1 and sigma_2
        are the electric conductivities of the two regions and
        grad(PotE_1) = (pot_1 - patch_value)/delta_1,
        grad(PotE_2) = (patch_value - pot_2)/delta_2 and
        UxB_1 and UxB_2 are the velocity contributions from the two regions.
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

        //- Default name of field on this region
        const word ePotName_ = "PotE";

        //- Default name of current density field
        const word Jname_ = "deltaJ";

        //- Default name of velocity contribution
        const word UxBname_ = "deltaUxB";

        //- Role name for distinguishing sides of an interface source
        const word terminalRole_;

        //bool coupled_ = false;


protected:

    //- Get useful values from patch
    virtual void getValues
    (
        tmp<scalarField>& sigma,
        tmp<scalarField>& sigmaByDelta,
        tmp<scalarField>& EPot
    ) const;

    //- Assign field to result which may have not been previously set
    void assign
    (
        tmp<scalarField>& result,
        const tmp<scalarField>& field
    ) const;

    //- Get suffix for complex cases (Re or Im)
    word suffix() const;

    //- Returns role of the terminal
    word getTerminalRole() const;

    //- Checks if source patch is valid
    bool isValidSource(word nbrRole) const;

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

    // Friend class
    //- Allow coupledCurrentDensityFvPatchVectorField to access private data
    //- and member functions to update current density from electric potential
    //friend class coupledCurrentDensityFvPatchVectorField;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
