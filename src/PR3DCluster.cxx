#include "WireCellData/PR3DCluster.h"

#include "TMatrixDEigen.h"
#include "TH2F.h"
#include "TVector3.h"


using namespace WireCell;

PR3DCluster::PR3DCluster(int cluster_id)
  : cluster_id(cluster_id)
{
  point_cloud = 0;
}

PR3DCluster::~PR3DCluster(){
  if (point_cloud!=(ToyPointCloud*)0)
    delete point_cloud;
}

// void AddCell(SlimMergeGeomCell* mcell, int *time_slices, int ntime_slice){
// }

void PR3DCluster::Create_point_cloud(){
  if (point_cloud!=(ToyPointCloud*)0)
    return;
  
  point_cloud = new ToyPointCloud();
  for (auto it = mcells.begin(); it!=mcells.end(); it++){
    SlimMergeGeomCell *mcell = (*it);
    PointVector pts = mcell->get_sampling_points();
    point_cloud->AddPoints(pts,mcell);
  }
  point_cloud->build_kdtree_index();
  //  std::cout << point_cloud->get_num_points() << std::endl;
  
}

std::pair<double,double> PR3DCluster::HoughTrans(Point&p , double dis){
  double theta, phi;
  TH2F *hough = new TH2F("","",180,0.,3.1415926,360,-3.1415926,3.1415926);
  double x0 = p.x;
  double y0 = p.y;
  double z0 = p.z;

  std::vector<std::pair<WireCell::SlimMergeGeomCell*,Point>>pts = point_cloud->get_closest_points(p,dis);
  double x,y,z,q;
  for (size_t i=0; i!=pts.size(); i++){
    x = pts.at(i).second.x;
    y = pts.at(i).second.y;
    z = pts.at(i).second.z;
    q = pts.at(i).first->get_q()/pts.at(i).first->get_sampling_points().size();
    if (q ==0) q = 1;
    TVector3 vec(x-x0,y-y0,z-z0);
    hough->Fill(vec.Theta(),vec.Phi(),q);
  }
  int maxbin = hough->GetMaximumBin();
  int a,b,c;
  hough->GetBinXYZ(maxbin,a,b,c);
  theta = hough->GetXaxis()->GetBinCenter(a);
  phi = hough->GetYaxis()->GetBinCenter(b);
  
  delete hough;
  return std::make_pair(theta,phi);
}


void PR3DCluster::Calc_PCA(){
  center.x=0; center.y=0; center.z=0;
  int nsum = 0;
  for (auto it = mcells.begin(); it!=mcells.end();it++){
    PointVector ps = (*it)->get_sampling_points();
    for (int k=0;k!=ps.size();k++){
      center.x += ps.at(k).x;
      center.y += ps.at(k).y;
      center.z += ps.at(k).z;
      nsum ++;
    }
  }

  for (int i=0;i!=3;i++){
    PCA_axis[i].x = 0;
    PCA_axis[i].y = 0;
    PCA_axis[i].z = 0;
  }
  
  if (nsum>=3){
    center.x /=nsum;
    center.y /=nsum;
    center.z /=nsum;
  }else{
    return;
  }
  TMatrixD cov_matrix(3,3);

  for (int i=0;i!=3;i++){
    for (int j=i;j!=3;j++){
      cov_matrix(i,j)=0;
      for (auto it = mcells.begin(); it!=mcells.end();it++){
	PointVector ps = (*it)->get_sampling_points();
	for (int k=0;k!=ps.size();k++){
	  if (i==0 && j==0){
	    cov_matrix(i,j) += (ps.at(k).x - center.x) * (ps.at(k).x - center.x);
	  }else if (i==0 && j==1){
	    cov_matrix(i,j) += (ps.at(k).x - center.x) * (ps.at(k).y - center.y);
	  }else if (i==0 && j==2){
	    cov_matrix(i,j) += (ps.at(k).x - center.x) * (ps.at(k).z - center.z);
	  }else if (i==1 && j==1){
	    cov_matrix(i,j) += (ps.at(k).y - center.y) * (ps.at(k).y - center.y);
	  }else if (i==1 && j==2){
	    cov_matrix(i,j) += (ps.at(k).y - center.y) * (ps.at(k).z - center.z);
	  }else if (i==2 && j==2){
	    cov_matrix(i,j) += (ps.at(k).z - center.z) * (ps.at(k).z - center.z);
	  }
	}
      }
    }
  }
  cov_matrix(1,0) = cov_matrix(0,1);
  cov_matrix(2,0) = cov_matrix(0,2);
  cov_matrix(2,1) = cov_matrix(1,2);
  
  TMatrixDEigen eigen(cov_matrix);
  TMatrixD eigen_values = eigen.GetEigenValues();
  TMatrixD eigen_vectors = eigen.GetEigenVectors();

  // std::cout << eigen_values(0,0) << " " << eigen_values(1,1) << " " << eigen_values(2,2) << std::endl;
  //std::cout << eigen_vectors(0,0) << " " << eigen_vectors(0,1) << " " << eigen_vectors(0,2) << std::endl;
  for (int i=0;i!=3;i++){
    PCA_axis[i].x = eigen_vectors(0,i)/sqrt(eigen_vectors(0,i)*eigen_vectors(0,i) + eigen_vectors(1,i)*eigen_vectors(1,i) + eigen_vectors(2,i)*eigen_vectors(2,i));
    PCA_axis[i].y = eigen_vectors(1,i)/sqrt(eigen_vectors(0,i)*eigen_vectors(0,i) + eigen_vectors(1,i)*eigen_vectors(1,i) + eigen_vectors(2,i)*eigen_vectors(2,i));
    PCA_axis[i].z = eigen_vectors(2,i)/sqrt(eigen_vectors(0,i)*eigen_vectors(0,i) + eigen_vectors(1,i)*eigen_vectors(1,i) + eigen_vectors(2,i)*eigen_vectors(2,i));

  }
  
  
  //std::cout << mean_x << " " << mean_y << " " << mean_z << std::endl;
}


void PR3DCluster::AddCell(SlimMergeGeomCell* mcell, int time_slice){
  if (cell_times_set_map.find(mcell)==cell_times_set_map.end()){
    std::set<int> times;
    times.insert(time_slice);
    cell_times_set_map[mcell]=times;
    mcells.push_back(mcell);
  }else{
    std::set<int>& times = cell_times_set_map[mcell];
    //if (find(times.begin(),times.end(),time_slice)==times.end()){
    times.insert(time_slice);
    // }
  }
  
  if (time_cells_set_map.find(time_slice)==time_cells_set_map.end()){
    SMGCSet mcells_1;
    mcells_1.insert(mcell);
    time_cells_set_map[time_slice] = mcells_1;
  }else{
    SMGCSet& mcells_1 = time_cells_set_map[time_slice];
    //if (find(mcells_1.begin(),mcells_1.end(), mcell) == mcells_1.end()){
    mcells_1.insert(mcell);
    //}
  }
}

void PR3DCluster::Remove_duplicated_mcells(){
  SMGCSelection to_be_removed;
  SMGCSelection to_be_saved;
  std::map<SlimMergeGeomCell*,SlimMergeGeomCell*> mcell_map;
  for (auto it = mcells.begin();it!=mcells.end();it++){
    SlimMergeGeomCell *mcell = (*it);
    bool save = true;
    for (auto it1=to_be_saved.begin(); it1!=to_be_saved.end(); it1++){
      SlimMergeGeomCell *mcell1 = *it1;
      if (mcell->IsSame(mcell1)){
	to_be_removed.push_back(mcell);
	mcell_map[mcell] = mcell1;
	save = false;
	break;
      }
    }
    if(save) to_be_saved.push_back(mcell);
  }
  // Now merge ...
  for (auto it = to_be_removed.begin(); it!= to_be_removed.end(); it++){
    SlimMergeGeomCell *mcell = *it; // to be removed
    SlimMergeGeomCell *mcell1 = mcell_map[mcell]; // repalced by this one

    // save time
    std::set<int> times = cell_times_set_map[mcell];
    for (auto it1 = times.begin(); it1!=times.end(); it1++){
      cell_times_set_map[mcell1].insert(*it1);
      time_cells_set_map[*it1].insert(mcell1);
      time_cells_set_map[*it1].erase(mcell);
    }

    // Now remove ...
    cell_times_set_map.erase(mcell);
    auto it2 = find(mcells.begin(), mcells.end(), mcell);
    mcells.erase(it2);
    delete mcell;
  }
  to_be_removed.clear();
  
  //std::cout << mcells.size() << " " << to_be_saved.size() << " " << to_be_removed.size() << std::endl;
}


SMGCSelection PR3DCluster::Is_Connected(PR3DCluster* cluster1 , int offset){
  SMGCSelection temp_mcells;

  for (auto it = mcells.begin(); it!= mcells.end();it++){
    SlimMergeGeomCell *mcell = *it;
    int time_slice = mcell->GetTimeSlice();
    std::map<int,SMGCSet>& time_mcell_map = cluster1->get_time_cells_set_map();
    // std::cout << time_slice << " " << time_mcell_map.size() << std::endl;
    // std::vector<int> times;
    // if (time_cells_set_map.find(time_slice)!=time_cells_set_map.end())
    //   times.push_back(time_slice);
    // if (time_cells_set_map.find(time_slice-1)!=time_cells_set_map.end())
    //   times.push_back(time_slice-1);
    // if (time_cells_set_map.find(time_slice+1)!=time_cells_set_map.end())
    //   times.push_back(time_slice+1);

    bool flag = false;
    // for (int i=0;i!=times.size();i++){
    if (time_mcell_map.find(time_slice)!=time_mcell_map.end()){
      SMGCSet& mcells1 = time_mcell_map[time_slice];
      for (auto it1 = mcells1.begin(); it1!=mcells1.end();it1++){
   	SlimMergeGeomCell *mcell1 = (SlimMergeGeomCell*)(*it1);
   	if (mcell->Overlap_fast(mcell1,offset)){
   	  flag = true;
   	  break;
   	}
      }
      //     //  if (flag) break;
    }
    if (flag) temp_mcells.push_back(mcell);
  }
  
  return temp_mcells;
}
