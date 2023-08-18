if __name__ == '__main__':
  
    from procare.open3d.open3d.registration import registration_icp
    from procare.open3d.open3d.registration import registration_ransac_based_on_feature_matching
    from procare.open3d.open3d.geometry import read_point_cloud
    from procare.open3d.open3d.registration import compute_cfpfh_feature
    from procare.open3d.open3d.geometry import estimate_normals
    from procare.open3d.open3d.geometry import KDTreeSearchParamHybrid
    from procare.open3d.open3d.registration import CorrespondenceCheckerBasedOnEdgeLength
    from procare.open3d.open3d.registration import CorrespondenceCheckerBasedOnDistance
    from procare.open3d.open3d.registration import RANSACConvergenceCriteria
    from procare.open3d.open3d.registration import ICPConvergenceCriteria
    from procare.open3d.open3d.registration import TransformationEstimationPointToPoint
    from procare.open3d.open3d.registration import TransformationEstimationPointToPlane

    from procare.open3d.open3d.visualization import draw_geometries

    from procare.convertpdb import _volsite_cavity_pdb_
    from procare.procarescores import _ph4_ext_

    import os
    import argparse
    import copy
    import numpy as np
    from time import strftime, localtime

    parser = argparse.ArgumentParser(description='Parameters for ProCare')

    # inputs
    parser.add_argument('-s', '--source', type=str,
        help='Source pdb file', 
        required=True)

    parser.add_argument('-t', '--target', type=str,
        help='Target pdb file', 
        required=True)

    args = parser.parse_args()


    FUNCTIONS = {
    'TransformationEstimationPointToPlane': TransformationEstimationPointToPlane,
    'TransformationEstimationPointToPoint': TransformationEstimationPointToPoint 
    #default: TransformationEstimationPointToPoint(with_scaling = False)
    }

    # Source and target loading
    """
    e.g.
    ##############################source_file##############################
    1ha3A.pcd
    ##############################source_prop##############################
    [[1, 'N'], [2, 'C'], [3, 'C'], [4, 'O'], [5, 'N'], [6, 'C'], [7, 'C'], [8, 'O'], [9, 'C'], [10, 'C']]
    ##############################source_color##############################
    [65280, 16711680, 16711680, 255, 65280, 16711680, 16711680, 255, 16711680, 16711680]

    """
    source_cavity = _volsite_cavity_pdb_()
    source_file, source_prop, source_color = source_cavity.pdb_to_pcd(args.source)
    if source_file != -1:
        # extract 3D coordinates
        source = read_point_cloud(source_file)
    
    target_cavity = _volsite_cavity_pdb_()
    target_file, target_prop, target_color = target_cavity.pdb_to_pcd(args.target)
    if target_file != -1:
        target = read_point_cloud(target_file)
    
    if source_file and target_file != -1:
        print("read in successfully!")




