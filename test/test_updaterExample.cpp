// use the catch file but do not add the #define CATCH_CONFIG_MAIN !!
#include "utility/extern/catchorg/catch2/catch.hpp"

#include <LeMonADE/core/Ingredients.h>
#include <LeMonADE/feature/FeatureMoleculesIO.h>

#include "updater/UpdaterExample.h"

typedef LOKI_TYPELIST_1(FeatureMoleculesIO) Features;
typedef ConfigureSystem<VectorInt3,Features,4> Config;
typedef Ingredients<Config> IngredientsType;

TEST_CASE( "run the updater" ) {
    // empty instance of ingredients
    IngredientsType ingredients;

    uint32_t myParameter(42);

    // just check for errors to be thrown
    UpdaterExample<IngredientsType> Hugo1 (ingredients, myParameter);
    REQUIRE_NOTHROW( Hugo1.initialize() );
    REQUIRE_NOTHROW( Hugo1.execute() );
    REQUIRE_NOTHROW( Hugo1.cleanup() );

    // check a condition (should return TRUE), for instance
    CHECK(myParameter == 42);

    // check something to return false
    CHECK_FALSE( myParameter == 2578472);

}