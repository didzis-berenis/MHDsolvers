/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2022 OpenFOAM Foundation
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
    Foam::harmonicElectromagneticModel

Description
    Base-class for fluid thermodynamic properties.

See also
    Foam::electromagneticModel

SourceFiles
    harmonicElectromagneticModel.C

\*---------------------------------------------------------------------------*/

#ifndef harmonicElectromagneticModel_H
#define harmonicElectromagneticModel_H

#include "electromagneticModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class harmonicElectromagneticModel Declaration
\*---------------------------------------------------------------------------*/

class harmonicElectromagneticModel
:
    virtual public electromagneticModel
{

    //- Electric potential [V]
        //- real
        volScalarField PotEre_;
        //- imaginary
        volScalarField PotEim_;

    //- Current density [A/m^2]
        //- real
        volVectorField Jre_;

        //- imaginary
        volVectorField Jim_;

    //- Magnetic field [T]
        //- real
        volVectorField Bre_;

        //- imaginary
        volVectorField Bim_;

    const word imaginaryPartName_ = "imaginary";

public:


    //- Runtime type information
    TypeName("harmonicElectromagneticModel");


    //- Declare run-time constructor selection table
    declareRunTimeNewSelectionTable
    (
        autoPtr,
        harmonicElectromagneticModel,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Constructors

        //- Construct from components
        harmonicElectromagneticModel
        (
            const fvMesh& mesh,
            const word& phaseName
        );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<harmonicElectromagneticModel> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );

        //- Disallow default bitwise copy construction
        harmonicElectromagneticModel(const harmonicElectromagneticModel&) = delete;

    //- Destructor
    virtual ~harmonicElectromagneticModel();


    // Member Functions

        // Access to thermodynamic state variables

        //- Electric potential [V]
        virtual volScalarField& PotE(const word& partName);

        //- Current density [A/m^2]
        virtual volVectorField& J(const word& partName);

        //- Magnetic field [T]
        virtual volVectorField& B(const word& partName);

        //- Add the given density correction to the current density field.
        virtual volVectorField& correctJ(const word& partName) = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Complied within electromagneticModels.C
#ifdef NoRepository
    #include "harmonicElectromagneticModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
