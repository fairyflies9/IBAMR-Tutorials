#ifndef included_update_target_point_positions
#define included_update_target_point_positions

#include <IBAMR_config.h>    //Added for IBAMR KD Module 2016-11 by NAB (1/27/17)
#include <IBTK_config.h>     //Added for IBAMR KD Module 2016-11 by NAB (1/27/17)
#include <SAMRAI_config.h>   //Added for IBAMR KD Module 2016-11 by NAB (1/27/17)

// Headers for basic PETSc functions
#include <petscsys.h>        //Added for IBAMR KD Module 2016-11 by NAB (1/27/17)

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>          //Added for IBAMR KD Module 2016-11 by NAB (1/27/17)
#include <CartesianGridGeometry.h>    //Added for IBAMR KD Module 2016-11 by NAB (1/27/17)
#include <LoadBalancer.h>             //Added for IBAMR KD Module 2016-11 by NAB (1/27/17)
#include <StandardTagAndInitialize.h> //Added for IBAMR KD Module 2016-11 by NAB (1/27/17)

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBExplicitHierarchyIntegrator.h>    //Added for IBAMR KD Module 2016-11 by NAB (1/27/17)
#include <ibamr/IBMethod.h>                         //Added for IBAMR KD Module 2016-11 by NAB (1/27/17)
#include <ibamr/IBStandardForceGen.h>               //Added for IBAMR KD Module 2016-11 by NAB (1/27/17)
#include <ibamr/IBStandardInitializer.h>            //Added for IBAMR KD Module 2016-11 by NAB (1/27/17)
#include <ibamr/INSCollocatedHierarchyIntegrator.h> //Added for IBAMR KD Module 2016-11 by NAB (1/27/17)
#include <ibamr/INSStaggeredHierarchyIntegrator.h>  //Added for IBAMR KD Module 2016-11 by NAB (1/27/17)
#include <ibtk/AppInitializer.h>                    //Added for IBAMR KD Module 2016-11 by NAB (1/27/17)
#include <ibtk/LData.h>                             //Added for IBAMR KD Module 2016-11 by NAB (1/27/17)
#include <ibtk/LDataManager.h>
#include <ibtk/muParserCartGridFunction.h>          //Added for IBAMR KD Module 2016-11 by NAB (1/27/17)
#include <ibtk/muParserRobinBcCoefs.h>              //Added for IBAMR KD Module 2016-11 by NAB (1/27/17)

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

/*
 * Update the positions of the target point specifications.
 */
void
update_target_point_positions(
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    LDataManager* l_data_manager,
    const double current_time,
    const double dt);

#endif //#ifndef included_update_target_point_positions
