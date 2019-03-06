/*--------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2013-2015 by 
  o/.|.\o    E   nvironment    | LeMonADE Principal Developers (see AUTHORS)
    ooo                        | 
----------------------------------------------------------------------------------

This file is part of LeMonADE.

LeMonADE is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

LeMonADE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with LeMonADE.  If not, see <http://www.gnu.org/licenses/>.

--------------------------------------------------------------------------------*/

// minimal Lemonade includes from lemonade library 
#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>
#include <LeMonADE/utility/TaskManager.h>
#include <LeMonADE/utility/RandomNumberGenerators.h>


// include your own updater
#include "updater/UpdaterExample.h"

// include your own analyze
#include "analyzer/AnalyzerExample.h"

// read in utilities
#include "utility/extern/catchorg/clara/clara.hpp"


int main(int argc, char* argv[])
{
    /* read arguments
    * +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
    */
    int32_t parameterAnalyzer(0);
    uint32_t parameterUpdater(0);

    bool showHelp;

    auto parser
    = clara::Opt( parameterAnalyzer, "parameterAnalyzer" )
        ["-a"]["--analyzer"]
        ("The parameter for the analyzer")
    | clara::Opt( parameterUpdater, "parameterUpdater" )
        ["-u"]["--updater"]
        ("The parameter for the updater")
    | clara::Help( showHelp );

    auto result = parser.parse( clara::Args( argc, argv ) );
    if( !result ) {
        std::cerr << "Error in command line: " << result.errorMessage() << std::endl;
        exit(1);
    }else if(showHelp == true){
        // add a description here
        std::cout << "lemonade example for Yves" << std::endl;

        // parameters are explained automatically
		parser.writeToStream(std::cout);
		exit(0);
    }
    /* initialize system
    * +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
    */
    // set the typelist
    typedef LOKI_TYPELIST_1(FeatureMoleculesIO) Features;
    // define maximal number of bonds
    const uint32_t max_bonds=4;
    // setup the molecules type
    typedef ConfigureSystem<VectorInt3,Features,max_bonds> Config;
    //setup ingredients type (3D Vektor, oben genannte features, maximale Anzahl Bindungspartner)
    typedef Ingredients<Config> IngredientsType;
    //create an instance of the Ingredients Type with the configuration above
    IngredientsType ingredients;

    /* set up random number generator (static object)
    * +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
    */
    RandomNumberGenerators rng;
    rng.seedAll();

    /* use the TaskManager
    * +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
    */
    TaskManager taskManager;
    // add your updater
    taskManager.addUpdater(new UpdaterExample<IngredientsType>(ingredients,parameterUpdater));
    // add your analyzer
    taskManager.addAnalyzer(new AnalyzerExample<IngredientsType>(ingredients,parameterAnalyzer));

    // run initialize and the execute function as often as your need, just once in this case
    taskManager.initialize();
    taskManager.run(1);
    taskManager.cleanup();

    /* */
    return true;
}
