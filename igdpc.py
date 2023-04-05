"""
The IGD^+-C indicator for preference-based evolutionary multi-objective optimization using a reference point
"""
import numpy as np
import click

def igdplus(ref_point_set, point_set):
    """   
    Calculate the IGD^+ value of given point sets.

    Parameters
    ----------
    ref_point_set: 2-d float numpy ndarray
        A reference point set for the IGD calculation
    point_set: 2-d float numpy ndarray
        A point set
    """        
    sum_dist = 0
    if point_set.ndim == 1:
        point_set = np.array([point_set])    
    
    for rp in ref_point_set:                    
        min_dist = np.inf
        for p in point_set:        
            tmp = 0
            for x, r in zip(p, rp):
                tmp += (max(x - r, 0))**2
            tmp_dist = np.sqrt(tmp)
            min_dist = min(min_dist, tmp_dist)                
        sum_dist += min_dist
    return sum_dist / len(ref_point_set)

def igd_c(igd_ref_point_set, point_set):
    """   
    Calculate the IGD-C value of given point sets.

    [Reference] A. Mohammadi, M. N. Omidvar, X. Li, and K. Deb, "Integrating user preferences and decomposition methods for many-objective optimization," in Proceedings of the IEEE Congress on Evolutionary Computation, CEC 2014, Beijing, China, July 6-11, 2014. IEEE, 2014, pp. 421-428. [Online]. Available: https://doi.org/10.1109/CEC.2014.6900595

    Parameters
    ----------
    pset_file_path_list: a list of file paths         
        A list of file paths to load the point sets 
    qi_file_path: a list of  file paths         
        A list of file paths to save the quality indicator values
    igd_refpset_file_path: a file path       
        A file path to load the reference point set for the IGD calculation
    pref_point: 1-d float numpy ndarray
        An "n_obj"-dimensional reference point
    roi_radius: float
        A radius of the ROI
    """    
    igd_value = igdplus(igd_ref_point_set, point_set)
    return igd_value
                           
@click.command()
@click.option('--igd_refpset_file_path', '-i', required=True, type=str, help='A file pah to the reference point set for IGD^+-C')
@click.option('--pset_file_path', '-p', required=True, type=str, help='A file pah to the reference point set for IGD^+-C')
def run(igd_refpset_file_path, pset_file_path):
    igd_ref_point_set = np.loadtxt(igd_refpset_file_path, delimiter=",", comments="#", dtype=float)
    if igd_ref_point_set.ndim == 1:
        igd_ref_point_set = np.array([igd_ref_point_set])    

    point_set = np.loadtxt(pset_file_path, delimiter=",", comments="#", dtype=float)
                
    score = igd_c(igd_ref_point_set, point_set)
    print("IGD^+-C(X) = {}".format(score))
    
if __name__ == '__main__':
    run()
