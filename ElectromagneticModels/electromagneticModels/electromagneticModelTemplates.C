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

\*---------------------------------------------------------------------------*/

#include "electromagneticModel.H"
#include "wordIOList.H"
#include "compileTemplate.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class ElectromagneticModel, class Table>
typename Table::iterator Foam::electromagneticModel::lookupCstrIter
(
    const dictionary& thermoDict,
    Table* tablePtr
)
{
    const word modelType(thermoDict.lookup("electromagneticType"));

    Info<< "Selecting electromagnetics model " << modelType << endl;

    typename Table::iterator cstrIter = tablePtr->find(modelType);

    if (cstrIter == tablePtr->end())
    {
        FatalErrorInFunction
            << "Unknown " << ElectromagneticModel::typeName << " type "
            << modelType << nl << nl
            << "Valid " << ElectromagneticModel::typeName << " types are:" << nl
            << tablePtr->sortedToc() << nl
            << exit(FatalError);
    }

    return cstrIter;
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class ElectromagneticModel>
Foam::autoPtr<ElectromagneticModel> Foam::electromagneticModel::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    const IOdictionary electromagneticDict
    (
        physicalProperties::findModelDict(mesh, phaseName)
    );

    typename ElectromagneticModel::fvMeshConstructorTable::iterator cstrIter =
        lookupCstrIter<ElectromagneticModel, typename ElectromagneticModel::fvMeshConstructorTable>
        (
            electromagneticDict,
            ElectromagneticModel::fvMeshConstructorTablePtr_
        );

    return autoPtr<ElectromagneticModel>(cstrIter()(mesh, phaseName));
}

// ************************************************************************* //
