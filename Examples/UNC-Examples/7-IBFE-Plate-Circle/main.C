// Copyright (c) 2002-2014, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for basic libMesh objects
#include <libmesh/boundary_info.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>
#include <ibtk/LData.h>                    //added for MODULE 2016-11 NAB
#include <ibtk/LDataManager.h>             //added for MODULE 2016-11 NAB

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

// Elasticity model data.
namespace ModelData
{
  // Problem parameters.
  static const double mu = 10.0; //Elastic modulus
  static double kappa_s = 1.0e6; // body force spring constant

 // void
 // coordinate1_mapping_function(
//			       libMesh::Point& X,
//			       const libMesh::Point& s,
//			       void* /*ctx*/)
 // {
 //   X(0) = s(1)+.2;
 //   X(1) = s(0);
 //   return;
 // }

  //void
  //coordinate2_mapping_function(
//			       libMesh::Point& X,
//			       const libMesh::Point& s,
//			       void* /*ctx*/)
//  {
//
//    X(0) = s(1)-.2;
//    X(1) = s(0);
//    return;
//  }

void
target1_force_function(
    VectorValue<double>& F,
	const TensorValue<double>& /*FF*/,
	const libMesh::Point& X,
	const libMesh::Point& s,
	Elem* const /*elem*/,
	const std::vector<const std::vector<double>*>& /*var_data*/,
	const std::vector<const std::vector<VectorValue<double> >*>& /*grad_var_data*/,
	double time,
    void* /*ctx*/)
  {

    libMesh::Point s_dump;
    if(time<1)      
      {
	s_dump(0)=s(0);
	s_dump(1)=s(1)+.25*time;
      }

    else
      {
	s_dump(1)=X(1);
	s_dump(0)=X(0);
      }
    F = kappa_s*(s_dump-X);
    return;
  }


  // Stress tensor functions.
void
PK1_dev_stress_function(
    TensorValue<double>& PP,
    const TensorValue<double>& FF,
	const libMesh::Point& /*X*/,
	const libMesh::Point& /*s*/,
	Elem* const /*elem*/,
	const std::vector<const std::vector<double>*>& /*var_data*/,
	const std::vector<const std::vector<VectorValue<double> >*>& /*grad_var_data*/,
	double /*time*/,
    void* /*ctx*/)
  {
    //const double r = sqrt((s(0) - 0.5)*(s(0) - 0.5));
    //if (r <= 0.1)
    //{
    //PP.zero();
    //}
    //    else
    //{
    PP = mu*FF;
    //}
    return;
  }// PK1_dev_stress_function

void
PK1_dil_stress_function(
    TensorValue<double>& PP,
    const TensorValue<double>& FF,
	const libMesh::Point& /*X*/,
	const libMesh::Point& /*s*/,
	Elem* const /*elem*/,
	const std::vector<const std::vector<double>*>& /*var_data*/,
	const std::vector<const std::vector<VectorValue<double> >*>& /*grad_var_data*/,
	double /*time*/,
    void* /*ctx*/)
  {
    //const double r = sqrt((s(0) - 0.5)*(s(0) - 0.5));
    //if (r <= 0.1)
    //{
    //  PP.zero();
    //}
    //else
    //  {
    PP = -mu*tensor_inverse_transpose(FF,NDIM);
    //  }    
    return;
  }// PK1_dil_stress_function
}
using namespace ModelData;

// Function prototypes
/*void
  output_data(
  Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
  Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
  Mesh& mesh,
  EquationSystems* equation_systems,
  const int iteration_num,
  const double loop_time,
  const string& data_dump_dirname);*/

/*******************************************************************************
 * For each run, the input filename and restart information (if needed) must   *
 * be given on the command line.  For non-restarted case, command line is:     *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 * For restarted run, command line is:                                         *
 *                                                                             *
 *    executable <input file name> <restart directory> <restart number>        *
 *                                                                             *
 *******************************************************************************/
int
main(
     int argc,
     char* argv[])
{
  // Initialize libMesh, PETSc, MPI, and SAMRAI.
  LibMeshInit init(argc, argv);
  SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
  SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
  SAMRAIManager::startup();

  {// cleanup dynamically allocated objects prior to shutdown

    // Parse command line options, set some standard options from the input
    // file, initialize the restart database (if this is a restarted run),
    // and enable file logging.
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
    Pointer<Database> input_db = app_initializer->getInputDatabase();

    // Get various standard options set in the input file.
    const bool dump_viz_data = app_initializer->dumpVizData();
    const int viz_dump_interval = app_initializer->getVizDumpInterval();
    const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();
    const bool uses_exodus = dump_viz_data && !app_initializer->getExodusIIFilename().empty();
	const string mesh1_exodus_filename = app_initializer->getExodusIIFilename("mesh1");
	const string mesh2_exodus_filename = app_initializer->getExodusIIFilename("mesh2");
    const bool dump_restart_data = app_initializer->dumpRestartData();
    const int restart_dump_interval = app_initializer->getRestartDumpInterval();
    const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

    const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
    const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
    const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
    if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
      {
	Utilities::recursiveMkdir(postproc_data_dump_dirname);
      }

    const bool dump_timer_data = app_initializer->dumpTimerData();
    const int timer_dump_interval = app_initializer->getTimerDumpInterval();

    // Create a simple FE mesh with Dirichlet boundary conditions.
    //
    // Note that boundary condition data must be registered with each FE
    // system before calling IBFEMethod::initializeFEData().
    //Mesh mesh1(NDIM);
    //Mesh mesh2(NDIM);
	Mesh mesh1(init.comm(), NDIM);
	Mesh mesh2(init.comm(), NDIM);

    const double dx = input_db->getDouble("DX");
    const double ds = input_db->getDouble("MFAC")*dx;
    string elem_type = input_db->getString("ELEM_TYPE");
	mesh1.read("beam.e");
	mesh1.prepare_for_use();
	mesh2.read("circle.e");
	mesh2.prepare_for_use();

	vector<MeshBase*> meshes(2);
    meshes[0] = &mesh1;
    meshes[1] = &mesh2;

    // Create major algorithm and data objects that comprise the
    // application.  These objects are configured from the input database
    // and, if this is a restarted run, from the restart database.
    Pointer<INSHierarchyIntegrator> navier_stokes_integrator;
    const string solver_type = app_initializer->getComponentDatabase("Main")->getString("solver_type");
    if (solver_type == "STAGGERED")
      {
	navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
								       "INSStaggeredHierarchyIntegrator", app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
      }
    else if (solver_type == "COLLOCATED")
      {
	navier_stokes_integrator = new INSCollocatedHierarchyIntegrator(
									"INSCollocatedHierarchyIntegrator", app_initializer->getComponentDatabase("INSCollocatedHierarchyIntegrator"));
      }
    else
      {
	TBOX_ERROR("Unsupported solver type: " << solver_type << "\n" <<
		   "Valid options are: COLLOCATED, STAGGERED");
      }

    Pointer<IBFEMethod> ib_method_ops = new IBFEMethod(
						       "IBFEMethod", app_initializer->getComponentDatabase("IBFEMethod"), meshes, app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"));
    Pointer<IBHierarchyIntegrator> time_integrator = new IBExplicitHierarchyIntegrator(
										       "IBHierarchyIntegrator", app_initializer->getComponentDatabase("IBHierarchyIntegrator"), ib_method_ops, navier_stokes_integrator);
    Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
											  "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>(
									      "PatchHierarchy", grid_geometry);
    Pointer<StandardTagAndInitialize<NDIM> > error_detector = new StandardTagAndInitialize<NDIM>(
												 "StandardTagAndInitialize", time_integrator, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
    Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
    Pointer<LoadBalancer<NDIM> > load_balancer = new LoadBalancer<NDIM>(
									"LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
    Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm = new GriddingAlgorithm<NDIM>(
										       "GriddingAlgorithm", app_initializer->getComponentDatabase("GriddingAlgorithm"), error_detector, box_generator, load_balancer);

    // Configure the IBFE solver.
	ib_method_ops->initializeFEEquationSystems();
    //FEDataManager* fe_data_manager = ib_method_ops->getFEDataManager();
    IBFEMethod::LagBodyForceFcnData target1_force_data(target1_force_function);
    IBFEMethod::PK1StressFcnData PK1_dev_stress_data(PK1_dev_stress_function);
    IBFEMethod::PK1StressFcnData PK1_dil_stress_data(PK1_dil_stress_function);
    PK1_dev_stress_data.quad_order = Utility::string_to_enum<libMeshEnums::Order>(input_db->getStringWithDefault("PK1_DEV_QUAD_ORDER","THIRD"));
    PK1_dil_stress_data.quad_order = Utility::string_to_enum<libMeshEnums::Order>(input_db->getStringWithDefault("PK1_DIL_QUAD_ORDER","FIRST"));
 
    //ib_method_ops->registerInitialCoordinateMappingFunction(coordinate1_mapping_function, 0);
    //ib_method_ops->registerInitialCoordinateMappingFunction(coordinate2_mapping_function, 1);
    ib_method_ops->registerLagBodyForceFunction(target1_force_data, 0);
    ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data, 0);
    ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data, 1);
    ib_method_ops->registerPK1StressFunction(PK1_dil_stress_data, 0);
    ib_method_ops->registerPK1StressFunction(PK1_dil_stress_data, 1);

        EquationSystems* mesh1_equation_systems = ib_method_ops->getFEDataManager(0)->getEquationSystems();
        EquationSystems* mesh2_equation_systems = ib_method_ops->getFEDataManager(1)->getEquationSystems();

    //EquationSystems* mesh1_equation_systems = ib_method_ops->getFEDataManager(0)->getEquationSystems();
    //EquationSystems* mesh2_equation_systems = ib_method_ops->getFEDataManager(1)->getEquationSystems();
 

    // Create Eulerian initial condition specification objects.
    if (input_db->keyExists("VelocityInitialConditions"))
      {
	Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
									"u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
	navier_stokes_integrator->registerVelocityInitialConditions(u_init);
      }

    if (input_db->keyExists("PressureInitialConditions"))
      {
	Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
									"p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
	navier_stokes_integrator->registerPressureInitialConditions(p_init);
      }

    // Create Eulerian boundary condition specification objects (when necessary).
    const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
    vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
    if (periodic_shift.min() > 0)
      {
	for (unsigned int d = 0; d < NDIM; ++d)
	  {
	    u_bc_coefs[d] = NULL;
	  }
      }
    else
      {
	for (unsigned int d = 0; d < NDIM; ++d)
	  {
	    ostringstream bc_coefs_name_stream;
	    bc_coefs_name_stream << "u_bc_coefs_" << d;
	    const string bc_coefs_name = bc_coefs_name_stream.str();

	    ostringstream bc_coefs_db_name_stream;
	    bc_coefs_db_name_stream << "VelocityBcCoefs_" << d;
	    const string bc_coefs_db_name = bc_coefs_db_name_stream.str();

	    u_bc_coefs[d] = new muParserRobinBcCoefs(
						     bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
	  }
	navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
      }

    // Create Eulerian body force function specification objects.
    if (input_db->keyExists("ForcingFunction"))
      {
	Pointer<CartGridFunction> f_fcn = new muParserCartGridFunction(
								       "f_fcn", app_initializer->getComponentDatabase("ForcingFunction"), grid_geometry);
	time_integrator->registerBodyForceFunction(f_fcn);
      }

    // Set up visualization plot file writers.
    Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
    if (uses_visit)
      {
	time_integrator->registerVisItDataWriter(visit_data_writer);
      }
	
	std::unique_ptr<ExodusII_IO> mesh1_exodus_io(uses_exodus ? new ExodusII_IO(mesh1) : NULL);
	std::unique_ptr<ExodusII_IO> mesh2_exodus_io(uses_exodus ? new ExodusII_IO(mesh2) : NULL);

    // Initialize hierarchy configuration and data on all patches.
    ib_method_ops->initializeFEData();
    time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

    // Deallocate initialization objects.
    app_initializer.setNull();

    // Print the input database contents to the log file.
    plog << "Input database:\n";
    input_db->printClassData(plog);

    // Write out initial visualization data.
    int iteration_num = time_integrator->getIntegratorStep();
    double loop_time = time_integrator->getIntegratorTime();
    if (dump_viz_data)
      {
	pout << "\n\nWriting visualization files...\n\n";
	if (uses_visit)
	  {
	    time_integrator->setupPlotData();
	    visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
	  }
	if (uses_exodus)
	  {
	    mesh1_exodus_io->write_timestep(mesh1_exodus_filename, *mesh1_equation_systems, iteration_num/viz_dump_interval+1, loop_time);
	    mesh2_exodus_io->write_timestep(mesh2_exodus_filename, *mesh2_equation_systems, iteration_num/viz_dump_interval+1, loop_time);
	  }
      }

    // Main time step loop.
    double loop_time_end = time_integrator->getEndTime();
    double dt = 0.0;
    while (!MathUtilities<double>::equalEps(loop_time,loop_time_end) &&
	   time_integrator->stepsRemaining())
      {
	iteration_num = time_integrator->getIntegratorStep();
	loop_time = time_integrator->getIntegratorTime();

	pout <<                                                    "\n";
	pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
	pout << "At beginning of timestep # " <<  iteration_num << "\n";
	pout << "Simulation time is " << loop_time              << "\n";

	dt = time_integrator->getMaximumTimeStepSize();
	time_integrator->advanceHierarchy(dt);
	loop_time += dt;

	pout <<                                                    "\n";
	pout << "At end       of timestep # " <<  iteration_num << "\n";
	pout << "Simulation time is " << loop_time              << "\n";
	pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
	pout <<                                                    "\n";

	// At specified intervals, write visualization and restart files,
	// print out timer data, and store hierarchy data for post
	// processing.
	iteration_num += 1;
	const bool last_step = !time_integrator->stepsRemaining();
	if (dump_viz_data && (iteration_num%viz_dump_interval == 0 || last_step))
	  {
	    pout << "\nWriting visualization files...\n\n";
	    if (uses_visit)
	      {
		time_integrator->setupPlotData();
		visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
	      }
	    if (uses_exodus)
	      {
		mesh1_exodus_io->write_timestep(mesh1_exodus_filename, *mesh1_equation_systems, iteration_num/viz_dump_interval+1, loop_time);
		mesh2_exodus_io->write_timestep(mesh2_exodus_filename, *mesh2_equation_systems, iteration_num/viz_dump_interval+1, loop_time);
                
	      }
	  }
	if (dump_restart_data && (iteration_num%restart_dump_interval == 0 || last_step))
	  {
	    pout << "\nWriting restart files...\n\n";
	    RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
	  }
	if (dump_timer_data && (iteration_num%timer_dump_interval == 0 || last_step))
	  {
	    pout << "\nWriting timer data...\n\n";
	    TimerManager::getManager()->print(plog);
	  }
	//            if (dump_postproc_data && (iteration_num%postproc_data_dump_interval == 0 || last_step))
	//{
	//    pout << "\nWriting state data...\n\n";
	//    output_data(patch_hierarchy,
	//                navier_stokes_integrator, mesh1, mesh1_equation_systems,
	//                iteration_num, loop_time, postproc_data_dump_dirname);
	//   output_data(patch_hierarchy,
	//                navier_stokes_integrator, mesh2, mesh2_equation_systems,
	//                iteration_num, loop_time, postproc_data_dump_dirname);
	//}
      }

    // Cleanup Eulerian boundary condition specification objects (when
    // necessary).
    for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

  }// cleanup dynamically allocated objects prior to shutdown
  
  SAMRAIManager::shutdown();
  return 0;
}// main

/*void
  output_data(
  Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
  Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
  Mesh& mesh,
  EquationSystems* equation_systems,
  const int iteration_num,
  const double loop_time,
  const string& data_dump_dirname)
  {
  plog << "writing hierarchy data at iteration " << iteration_num << " to disk" << endl;
  plog << "simulation time is " << loop_time << endl;

  // Write Cartesian data.
  string file_name = data_dump_dirname + "/" + "hier_data.";
  char temp_buf[128];
  sprintf(temp_buf, "%05d.samrai.%05d", iteration_num, SAMRAI_MPI::getRank());
  file_name += temp_buf;
  Pointer<HDFDatabase> hier_db = new HDFDatabase("hier_db");
  hier_db->create(file_name);
  VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
  ComponentSelector hier_data;
  hier_data.setFlag(var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getVelocityVariable(), navier_stokes_integrator->getCurrentContext()));
  hier_data.setFlag(var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getPressureVariable(), navier_stokes_integrator->getCurrentContext()));
  patch_hierarchy->putToDatabase(hier_db->putDatabase("PatchHierarchy"), hier_data);
  hier_db->putDouble("loop_time", loop_time);
  hier_db->putInteger("iteration_num", iteration_num);
  hier_db->close();

  // Write Lagrangian data.
  file_name = data_dump_dirname + "/" + "fe_mesh.";
  sprintf(temp_buf, "%05d", iteration_num);
  file_name += temp_buf;
  file_name += ".xda";
  mesh.write(file_name);
  file_name = data_dump_dirname + "/" + "fe_equation_systems.";
  sprintf(temp_buf, "%05d", iteration_num);
  file_name += temp_buf;
  equation_systems->write(file_name, (EquationSystems::WRITE_DATA | EquationSystems::WRITE_ADDITIONAL_DATA));
  return;
  }// output_data
*/
