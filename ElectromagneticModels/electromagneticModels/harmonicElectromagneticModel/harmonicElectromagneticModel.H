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
        const char* JreName_ = "Jre";

        //- imaginary
        const char* JimName_ = "Jim";

    //- Magnetic field [T]
        //- real
        const char* BreName_ = "Bre";

        //- imaginary
        const char* BimName_ = "Bim";

public:


    //- Runtime type information
    TypeName("harmonicElectromagneticModel");


    //- Declare run-time constructor selection table
/*
    declareRunTimeNewSelectionTable
    (
        autoPtr,
        harmonicElectromagneticModel,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );
*/

    // Constructors

        //- Construct from components
        harmonicElectromagneticModel
        (
            const fvMesh& mesh,
            const word& phaseName
        );


    // Selectors
/*
        //- Standard selection based on fvMesh
        static autoPtr<harmonicElectromagneticModel> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );
*/
        //- Disallow default bitwise copy construction
        harmonicElectromagneticModel(const harmonicElectromagneticModel&) = delete;

    //- Destructor
    virtual ~harmonicElectromagneticModel();


    // Member Functions

        // Access to thermodynamic state variables

        //- Electric potential [V]
        virtual const volScalarField& PotE(bool imaginary=false);

        //- Current density [A/m^2]
        virtual volVectorField& J(bool imaginary=false);

        //- Magnetic field [T]
        virtual volVectorField& B(bool imaginary=false);

        //- Add the given density correction to the current density field.
        //virtual volVectorField& correctJ(bool imaginary=false) = 0;
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
