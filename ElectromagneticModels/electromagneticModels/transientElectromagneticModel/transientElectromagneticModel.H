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
    Foam::transientElectromagneticModel

Description
    Base-class for fluid thermodynamic properties.

See also
    Foam::electromagneticModel

SourceFiles
    transientElectromagneticModel.C

\*---------------------------------------------------------------------------*/

#ifndef transientElectromagneticModel_H
#define transientElectromagneticModel_H

#include "electromagneticModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class transientElectromagneticModel Declaration
\*---------------------------------------------------------------------------*/

class transientElectromagneticModel
:
    virtual public electromagneticModel
{

    //- Pressure [Pa]
    volScalarField PotE_;

    //- Compressibility [s^2/m^2]
    volVectorField J_;

    //- Compressibility [s^2/m^2]
    volVectorField B_;
public:


    //- Runtime type information
    TypeName("transientElectromagneticModel");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        transientElectromagneticModel,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Selectors

        //- Standard selection based on fvMesh
        static autoPtr<transientElectromagneticModel> New
        (
            const fvMesh&,
            const word& phaseName=word::null
        );


    //- Destructor
    virtual ~transientElectromagneticModel();


    // Member Functions

        // Access to thermodynamic state variables

            //- Pressure [Pa]
            //  Non-const access allowed for transport equations
            //virtual volScalarField& p() = 0;

            //- Pressure [Pa]
            virtual volScalarField& PotE() const = 0;

            //- Compressibility [s^2/m^2]
            virtual volVectorField& J() const = 0;

            //- Compressibility [s^2/m^2]
            virtual volVectorField& B() const = 0;

            //- Add the given density correction to the density field.
            //  Used to update the density field following pressure solution
            virtual volVectorField& correctJ() = 0;


        // Fields derived from thermodynamic state variables
/*
            //- Rename the thermodynamic density field if stored and return
            //  This is used by solvers which create a separate continuity rho
            //  [kg/m^3]
            virtual tmp<volScalarField> renameRho() = 0;

            //- Gamma = Cp/Cv []
            virtual tmp<volScalarField> gamma() const = 0;

            //- Gamma = Cp/Cv for patch []
            virtual tmp<scalarField> gamma
            (
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Molecular weight [kg/kmol]
            virtual tmp<volScalarField> W() const = 0;

            //- Molecular weight for patch [kg/kmol]
            virtual tmp<scalarField> W(const label patchi) const = 0;


        // Access to transport state variables

            //- Dynamic viscosity of mixture [kg/m/s]
            virtual tmp<volScalarField> mu() const = 0;

            //- Dynamic viscosity of mixture for patch [kg/m/s]
            virtual tmp<scalarField> mu(const label patchi) const = 0;

            //- Kinematic viscosity of mixture [m^2/s]
            virtual tmp<volScalarField> nu() const;

            //- Kinematic viscosity of mixture for patch [m^2/s]
            virtual tmp<scalarField> nu(const label patchi) const;
*/
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
