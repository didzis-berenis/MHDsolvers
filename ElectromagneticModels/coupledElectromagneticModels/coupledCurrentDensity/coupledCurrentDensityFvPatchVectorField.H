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
    Mixed boundary condition for temperature, to be used for heat-transfer
    with another region in a CHT case.  Optional thin wall material layer
    resistances can be specified through thicknessLayers and kappaLayers
    entries.

    Specifies gradient and temperature such that the equations are the same
    on both sides:
      - refGradient = ePotByDeltaNbr*sigmaNbr/sigma
      - refValue = neighbour value
      - valueFraction = 0, but switches to 1 for small sigma

    where ePotByDeltaNbr is the gradient
    and sigmaNbr is the conductivity of neighbour patch.

Usage
    \table
        Property     | Description             | Required    | Default value
        ePot         | name of the field       | no          | ePot
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            coupledTemperature;
    }
    \endverbatim

    Needs to be on underlying mapped(Wall)FvPatch.

    The patch thermal conductivity \c kappa is obtained from the region
    thermophysicalTransportModel so that this boundary condition can be applied
    directly to either fluid or solid regions.

    Note that in order to provide an optional heat source either qs or Qs
    should be specified, not both.

SourceFiles
    coupledCurrentDensityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef coupledCurrentDensityFvPatchVectorField_H
#define coupledCurrentDensityFvPatchVectorField_H

#include "directionMixedFvPatchFields.H"

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

        bool coupled_ = false;


protected:

    //- Get the patch sigma
    virtual void getThis
    (
        tmp<scalarField>& sigma
    ) const;

    //- Add field to result which may have not been previously set
    void add
    (
        tmp<scalarField>& result,
        const tmp<scalarField>& field
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
