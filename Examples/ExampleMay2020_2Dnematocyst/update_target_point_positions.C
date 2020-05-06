#include "update_target_point_positions.h"
#include <ibamr/IBTargetPointForceSpec.h>

void
update_target_point_positions(
			      tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
			      LDataManager* const l_data_manager,
			      const double current_time,
			      const double dt)
{
  const int finest_ln = hierarchy->getFinestLevelNumber();

  static const double pi = 4*atan(1);
  static const double V = 18; //velocity of the wing during translation (meters/sec)
  static const double time_to_acc = 0.1*0.000000722*5/1.3; //time to accelerate up to V (s)

  static const double L1 = 0.0005; // length of computational domain (meters)
  static const int N1 = 512; // number of cartesian grid meshwidths at the finest level of the AMR grid

  // Find out the Lagrangian index ranges.
  const std::pair<int,int>& plate2d_lag_idxs = l_data_manager->getLagrangianStructureIndexRange(0, finest_ln);
  const std::pair<int,int>& rubberband2d_lag_idxs = l_data_manager->getLagrangianStructureIndexRange(1, finest_ln);
	
  // Get the LMesh (which we assume to be associated with the finest level of
  // the patch hierarchy).  Note that we currently need to update both "local"
  // and "ghost" node data.
  Pointer<LMesh> mesh = l_data_manager->getLMesh(finest_ln);
  vector<LNode*> nodes;
  nodes.insert(nodes.end(), mesh->getLocalNodes().begin(), mesh->getLocalNodes().end());
  nodes.insert(nodes.end(), mesh->getGhostNodes().begin(), mesh->getGhostNodes().end());

  // Update the target point positions in their associated target point force
  // specs.
  tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(finest_ln);
  for (vector<LNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
    {
      LNode* node_idx = *it;
      IBTargetPointForceSpec* force_spec = node_idx->getNodeDataItem<IBTargetPointForceSpec>();
      if (force_spec == NULL) continue;  // skip to next node

      // Here we update the position of the target point.
      //
      // NOTES: lag_idx      is the "index" of the Lagrangian point (lag_idx = 0, 1, ..., N-1, where N is the number of Lagrangian points)
      //        X_target     is the target position of the target point
      //        X_target[0]  is the x component of the target position
      //        X_target[1]  is the y component of the target position
      //        X_target[2]  is the z component of the target position (for a 3D simulation)
      //
      // The target position is shifted to the left or right by the
      // increment dt*V

      const int lag_idx = node_idx->getLagrangianIndex();
      //Depending on the version of IBAMR, you need to select one of the ways of accessing target point positions. 
      //TinyVector<double,NDIM>& X_target = force_spec->getTargetPointPosition();
      //IBTK::Vector<double,NDIM>& X_target = force_spec->getTargetPointPosition();
      Point& X_target = force_spec->getTargetPointPosition();



      if (current_time <= time_to_acc)
	{
	  if (plate2d_lag_idxs.first <= lag_idx && lag_idx < plate2d_lag_idxs.second)
	    {
	      X_target[1]+=current_time/time_to_acc*V*dt;
	    }
	}
      else 
	{
	  if (plate2d_lag_idxs.first <= lag_idx && lag_idx < plate2d_lag_idxs.second)
	    {
	      X_target[1]+=V*dt;
	    }
	}
    }
  return;
}// update_target_point_positions
