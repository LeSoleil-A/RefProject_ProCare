# ----------------------------------------------------------------------------
# <                     ProCare For PDB only use ICP                         >
#                            Seems not effective
# ----------------------------------------------------------------------------

# Step One Function
def process_pointcloud(pointcloud_, radius_normal_, 
        radius_feature_, max_nn_normal_, max_nn_feature_):

    """This fuction estimates normals and calculate 
    the fast point feature histogramm"""

    # returns bool
    estimate_normals(pointcloud_, KDTreeSearchParamHybrid(
               radius=radius_normal_, max_nn=max_nn_normal_)) 
    # compute FPFH feature for a point cloud
    cfpfh = compute_cfpfh_feature(pointcloud_, KDTreeSearchParamHybrid(
               radius=radius_feature_, max_nn=max_nn_feature_))
    return cfpfh, pointcloud_

# Step Two Function
def fine_registration(source_, target_, distance_threshold_, 
        transformation_type_, relative_rmse_, relative_fitness_, max_iter_):
    
    function_transtype = FUNCTIONS[transformation_type_]
    # default TransformationEstimationPointToPoint: with_scaling = False
    # Function for ICP registration
    result = registration_icp(source_, target_,
        max_correspondence_distance=distance_threshold_,
        estimation_method=function_transtype(),
        criteria=ICPConvergenceCriteria(relative_fitness_, relative_rmse_, max_iter_))
    return result

# Version 1.0: Step Three Function
def transform(pdb_ofile_, transformed_coords_, source_color_):
    rotated_pdb = _volsite_cavity_pdb_()
    rot_coords = []
    for point, color in zip(transformed_coords_, source_color_):
        rot_coords.append([point[0], point[1], point[2], color])

    rotated_pdb.write_pdb_fake(pdb_ofile_, rot_coords)

# Version 2.0: Step Three Function
def transform_coor(source_file_, pdb_ofile_, transformed_coords_, source_color_):
    rotated_pdb = _volsite_cavity_pdb_()
    rot_coords = []
    for point, color in zip(transformed_coords_, source_color_):
        rot_coords.append([point[0], point[1], point[2], color])

    rotated_pdb.write_pdb_from_coor(source_file_, pdb_ofile_, rot_coords)



if __name__ == '__main__':
  
    from procare.open3d.open3d.registration import registration_icp
    from procare.open3d.open3d.geometry import read_point_cloud
    from procare.open3d.open3d.registration import compute_cfpfh_feature
    from procare.open3d.open3d.geometry import estimate_normals
    from procare.open3d.open3d.geometry import KDTreeSearchParamHybrid
    from procare.open3d.open3d.registration import ICPConvergenceCriteria
    from procare.open3d.open3d.registration import TransformationEstimationPointToPoint
    from procare.open3d.open3d.registration import TransformationEstimationPointToPlane

    from procare.open3d.open3d.visualization import draw_geometries

    from procare.convertpdb import _volsite_cavity_pdb_
    from procare.procarescorespdb import _ph4_ext_pdb_

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

    # normals
    parser.add_argument('-nr', '--normalrad', type=float,
        help='Radius for local surface normal estimation on a point', 
        required=False,
        default=3.1)

    parser.add_argument('-nm', '--normalmaxn', type=int,
        help=('Maximum number of neighbors to consider for local surface normal '
              'estimation on a point'), 
        required=False,
        default=471)

    # features
    parser.add_argument('-fr', '--featurerad', type=float,
        help='Radius for local surface feature estimation on a point', 
        required=False,
        default=3.1)

    parser.add_argument('-fm', '--featuremaxn', type=int,
        help=('Maximum number of neighbors to consider for local surface feature '
              'estimation on a point'), 
        required=False,
        default=135)       

    # icp
    parser.add_argument('-it', '--icptranstype', type=str,
        help='Transformation estimation type for ICP registration', 
        required=False,
        default='TransformationEstimationPointToPoint')

    parser.add_argument('-id', '--icpdist', type=float,
        help='Distance threshold for correspondences set in ICP registration', 
        required=False,
        default=3)

    parser.add_argument('-ir', '--icprmse', type=float,
        help='RMSE relative threshold for ICP terminaison', # default observed as acceptable
        required=False,
        default=10e-6)

    parser.add_argument('-if', '--icpfitness',type=float,
        help='Fitness relative threshold for ICP terminaison', 
        required=False,
        default=10e-6)

    parser.add_argument('-ii', '--icpiter', type=int,
        help='ICP convergence criteria: maximum iteration', 
        required=False,
        default=100)     
    
    # transform output
    parser.add_argument('--transform', action='store_true',
        help='output rotated mol2', 
        required=False)
    
    parser.add_argument('--ligandtransform', type=str, nargs='+',
        help='output rotated ligand and/or protein mol2', 
        required=False)
    
    # similarity output
    parser.add_argument('-so', '--scoreoutput', type=str,
        help='Simplifiled output file: scores', 
        required=False,
        default='procare_scores.tsv')
    
    parser.add_argument('-o', '--output', type=str,
        help='Complete output file', 
        required=False,
        default='procare.tsv')
    
    parser.add_argument('-p', '--paramid', type=str,
        help='ID for parameters identification', 
        required=False,
        default='default')

    parser.add_argument('-c', '--classification', type=str,
        help='Class for retrospective screening: 0 or 1',
        choices=['0', '1'],
        required=False,
        default='NAN')



    args = parser.parse_args()


    FUNCTIONS = {
    'TransformationEstimationPointToPlane': TransformationEstimationPointToPlane,
    'TransformationEstimationPointToPoint': TransformationEstimationPointToPoint 
    #default: TransformationEstimationPointToPoint(with_scaling = False)
    }

    # Source and target loading
    # sequential loading: sourcefile and target files can have the same file names
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
        # print("read in successfully!")

        ########################################### STEP 1: Read in the point clouds#########################################

        # estimates normals and calculate the fast point feature histogramm (FPFH)
        # source_cfpfh: FPFH feature for a point cloud (open3d.pipelines.registration.Feature)
        # source: point cloud (with normals?)

        source_cfpfh, source = process_pointcloud(pointcloud_=source,
                                             radius_normal_=args.normalrad,
                                             radius_feature_=args.featurerad,
                                             max_nn_normal_=args.normalmaxn,
                                             max_nn_feature_=args.featuremaxn)
        
        target_cfpfh, target = process_pointcloud(pointcloud_=target,
                                         radius_normal_=args.normalrad,
                                         radius_feature_=args.featurerad,
                                         max_nn_normal_=args.normalmaxn,
                                         max_nn_feature_=args.featuremaxn)
        
        # if source != None:
        #     print("process success")

        ########################################### STEP 2: ICP alignement ##########################################
        # Initial ICP alignement based of features
        """
        Input:
            --icprmse: RMSE relative threshold for ICP terminaison
            --icpfitness: Fitness relative threshold for ICP terminaison
        Output:
            result_fine_cfpfh: contains the registration results (open3d.pipelines.registration.RegistrationResult)
        """         
        result_fine_cfpfh = fine_registration(source_=source,
                                        target_=target,
                                        distance_threshold_=args.icpdist,
                                        transformation_type_=args.icptranstype,
                                        relative_rmse_=args.icprmse,
                                        relative_fitness_=args.icpfitness,
                                        max_iter_=args.icpiter)     
        
        # if result_fine_cfpfh != None:
        #     print("fine registration success")
        print("result_fine_cfpfh: ", result_fine_cfpfh.transformation)

        ########################################### STEP 3: Calculate Similarity##########################################
        source_transformed_cfpfh = copy.deepcopy(source)

        # open3d.geometry.Geometry3D (The estimated transformation matrix: 4*4 float64 numpy array)
        """
        Is this what we need?
        e.g:result_fine_cfpfh.transformation: 
            [[ 9.99820290e-01  1.87492883e-02  2.80212393e-03  5.28324300e+01]
             [-1.87562194e-02  9.99821041e-01  2.46803150e-03 -9.68504437e-01]
             [-2.75534863e-03 -2.52014522e-03  9.99993028e-01  1.74062308e-01]
             [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  1.00000000e+00]]
        """
        # source_transformed_cfpfh: transformation result
        source_transformed_cfpfh.transform(result_fine_cfpfh.transformation)

        if args.transform:
            rot_file_rot = 'rot_{}.pdb'.format(os.path.splitext(source_file)[0])
            transform_coor(source_file_=args.source,
                      pdb_ofile_=rot_file_rot,
                      transformed_coords_=source_transformed_cfpfh.points,
                      source_color_=source_color)              

        ph4_ext = _ph4_ext_pdb_(source_transformed_cfpfh.points, target.points, 
                                            source_prop, target_prop, 1.5)
        
        # print(np.asarray(target.points))

        ratio_aligned, \
            ratio_C_in_aligned, ratio_N_in_aligned, \
            ratio_O_in_aligned, ratio_S_in_aligned = ph4_ext.get_similarity_by_rules()
        
        # Use tversky_similarity to calculate similarity score
        score = ph4_ext.tversky_similarity()


        ########################################### Write out the result ##########################################

        if not os.path.isfile(args.scoreoutput):
            with open(args.scoreoutput, "w") as of:
                of.write("Source\tTarget\tScore\ns")
        
        with open(args.scoreoutput, "a") as of:
            of.write("{}\t{}\t{}\n".format(os.path.splitext(source_file)[0],
                                           os.path.splitext(target_file)[0],
                                           score))
            
        # output contributions of the differnt ph4 to the global score
        # output matrix components

        if not os.path.isfile(args.output):
            with open(args.output, "w") as of:
                of.write("Param_id\tSource\tTarget\tClass\tScore\t"

                         "C_contrib\tN_contrib\tO_contrib\tS_contrib\t"

                         "G_fitness\tG_RMSE\tICP_fitness\tICP_RMSE\t"

                         "ICP_11\tICP_12\tICP_13\tICP_14\t"
                         "ICP_21\tICP_22\tICP_23\tICP_24\t"
                         "ICP_31\tICP_32\tICP_33\tICP_34\t"
                         "ICP_41\tICP_42\tICP_43\tICP_44\n")    

        with open(args.output, 'a') as of:
            of.write(("{}\t{}\t{}\t{}\t{}\t"
                      "{}\t{}\t{}\t{}\t"
                      "{}\t{}\t{}\t{}\t"
                      "{}\t{}\t{}\t{}\t"
                      "{}\t{}\t{}\t{}\t"
                      "{}\t{}\t{}\t{}\t"
                      "{}\t{}\t{}\t{}\t"
                      "{}\t{}\t{}\t{}\t"
                      "{}\t{}\t{}\t{}\t"
                      "{}\t{}\t{}\t{}\t"
                      "{}\t{}\t{}\t{}\n"
                      ).format(args.paramid,
                               os.path.splitext(source_file)[0],
                               os.path.splitext(target_file)[0],
                               args.classification,
                               score,
                               ratio_C_in_aligned,
                               ratio_N_in_aligned,
                               ratio_O_in_aligned,
                               ratio_S_in_aligned,
                               result_fine_cfpfh.fitness,
                               result_fine_cfpfh.inlier_rmse,
                               np.array(result_fine_cfpfh.transformation)[0][0],
                               np.array(result_fine_cfpfh.transformation)[0][1],
                               np.array(result_fine_cfpfh.transformation)[0][2],
                               np.array(result_fine_cfpfh.transformation)[0][3],
                               np.array(result_fine_cfpfh.transformation)[1][0],
                               np.array(result_fine_cfpfh.transformation)[1][1],
                               np.array(result_fine_cfpfh.transformation)[1][2],
                               np.array(result_fine_cfpfh.transformation)[1][3],
                               np.array(result_fine_cfpfh.transformation)[2][0],
                               np.array(result_fine_cfpfh.transformation)[2][1],
                               np.array(result_fine_cfpfh.transformation)[2][2],
                               np.array(result_fine_cfpfh.transformation)[2][3],
                               np.array(result_fine_cfpfh.transformation)[3][0],
                               np.array(result_fine_cfpfh.transformation)[3][1],
                               np.array(result_fine_cfpfh.transformation)[3][2],
                               np.array(result_fine_cfpfh.transformation)[3][3]
                               ))

        # # cleaning
        # # if source_file == target_file:
        # #     os.system('rm {}'.format(source_file))
        # # else:
        # #     os.system('rm {} {}'.format(source_file, target_file))
