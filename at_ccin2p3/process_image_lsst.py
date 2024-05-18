import os
import shlex, subprocess
from astropy.io import fits
import file_conf as fc
import glob
import time
import psutil
import numpy as np
import stack_for_stage1 as stack_lib
import butler4euclid as butler
import eas

"""
ProcessStage1LSST is the main class to processing LSST data for Euclid stage 1
"""

def analyse_name_ccd_file(p_file):
    '''Extract visit, raft and sensor from name file raw image
    analyze name raw image LSST format : lsst_e_601482_f2_R33_S00_E000.fits
    :param p_file:
    '''
    out_dict = {}
    name = p_file.split('_')
    out_dict['visit'] = name[2]
    out_dict['raft'] = name[4][1] + "," + name[4][2]
    out_dict['sensor'] = name[5][1] + "," + name[5][2]
    return out_dict

    
def extract_visit_from_xml_name(xml_name):
    '''
    format of xml name used
    EUC_EXT_DPDEXTSINGLEEPOCHFRAME_LSST-370990-R33-S01_20210601T081226.348Z        
    '''
    return xml_name.split('-')[1]


def keep_name(a):
    '''
    return 'name' if a='toto/tutu/name.xml'
    :param a:
    '''    
    return a.split('/')[-1].split('.xml')[0]


def div_list(l_file, nb_proc):
    '''divide list in nb_proc chunck
    example return: 
      [[1, 2, 3, 4], [5, 6, 7], [8, 9, 10]]
    for
      l_file = [1,2,3,4,5,6,7,8,9,10]
      nb_proc = 3
      
    :param l_file:
    :param nb_proc:
    '''
    n_file = len(l_file)
    nb_file_by_chunck = n_file // nb_proc
    # all chunk has nb_file_by_chunck elements
    idx = np.ones(nb_proc + 1, dtype=np.int32) * nb_file_by_chunck
    # rest on the first chunk: => add 1  
    idx[:n_file % nb_proc + 1] += 1
    idx[0] = 0
    # create chunk index in l_file
    idx = np.cumsum(idx)
    # print(idx)
    l_list = [l_file[a:b] for a, b in zip(idx[:-1], idx[1:])]
    # print(l_list)
    return l_list


def write_file_from_list(m_list, name_file):
    with open(name_file, 'w') as d_file:
        for elt in m_list:
            d_file.write(elt + "\n")
                                                        

class ManageSetButlers(object):
    
    def __init__(self, root_butler):
        self.root = root_butler
        
    def _keep_only_file_no_extension(self, l_xml):
        l_only_xml = list(map(keep_name, l_xml))
        return l_only_xml
            
    def get_all_xml(self):
        cur_dir = os.getcwd()
        os.chdir(self.root)
        l_xml = glob.glob('butler_*/stage1/*.xml')
        os.chdir(cur_dir)
        l_xml = self._keep_only_file_no_extension(l_xml)
        return l_xml
    
    def _sef_for_eas(self, eas_obj):
        assert isinstance(eas_obj, eas.EasToolsToIngest)
        in_eas = eas_obj.get_list_sef_in_eas()
        in_but = self.get_all_xml()
        to_eas = list(set(in_but) - set(in_eas))
        print("in eas: ", len(in_eas))
        print("in cc : ", len(in_but))        
        self.sef4eas = to_eas
        self.nb_sef_in_eas = len(in_eas)
        self.nb_sef_in_but = len(in_but)
        del in_eas
        del in_but
        
    def load_balancing_eas(self, nb_proc, eas_obj):
        # create lb_eas directory
        lb_dir = os.path.join(self.root, "lb_eas")
        try:
            os.mkdir(lb_dir)
        except FileExistsError:
            pass
        # clean old load balancing
        cmd = f"rm -rf {lb_dir}/lb_*.txt"
        os.system(cmd)        
        # create nb_proc file with xml to ingest
        self._sef_for_eas(eas_obj)
        self.nb_sef_job = 0      
        if len(self.sef4eas) == 0:
            if self.nb_sef_in_but == 0:
                print("No SEF xml in butler_xx/stage1 directory !? === > Create it before call this script")
            else:
                print('Seems all SEF are in EAS')
        else:
            # divide SEF in nb_proc list
            lb_eas = div_list(self.sef4eas, nb_proc)
            self.nb_sef_job = len(lb_eas[0])
            print(f"{len(self.sef4eas)} SEF to transfert to EAS.")
            print(f"{nb_proc} jobs to ingest ~{self.nb_sef_job} SEF each.")
            for idx, sef_set in enumerate(lb_eas):
                file_sef = os.path.join(lb_dir, f"lb_eas_{idx}.txt")          
                with open(file_sef, 'w') as f_sef:          
                    for sef in sef_set:
                        f_sef.write(sef + "\n")
                        
class ProcessStage1LSST(object):
    """Process image with same id visit or reference star catalog used for simulation
    """
    
    def __init__(self):        
        self.rank = 0
        self.nb_proc = 10        
        self.dir_lb = 'lb_image'
        self.file_lb = 'list_image_for_rank_xx.txt'
        self.dir_butler = 'butler_x'
        self.limit_max_file = 10000
        self.skip_proc = False
    
    def init_with_file_param(self, p_params):
        assert isinstance(p_params, fc.FileConfLSSTProcessing)
        print(p_params)
        self.file_name_image = p_params.file_list_images        
        self.skip_proc = p_params.skip_proc
        self.params = p_params
        self.stack = stack_lib.factory_for_stage1_stack(self.params.singularity_lsst_stack)
                
    def init_env_var(self):        
        self.rank = int(os.environ["SGE_TASK_ID"]) - 1
        # self.nb_proc = int(os.environ["SGE_TASK_LAST"])        
        self.tmp_dir = os.environ["TMPDIR"]
        self.file_name_image = os.path.join(self.params.root_butlers, self.dir_lb, self.file_lb.replace('xx', str(self.rank)))        
    
    def extract_visit_ids_from_name(self, fname):
        """ 
        return id visitstring '100154' from raw image name
        
        Example:
            >>> fname="SIM_SC8_GSIR_SWF1_RUBIN_PPO_7/data/EUC_SIM_LSST-R14S00-1492656-SCIENCE-z_07CD30AD0911-0042104_20201021T162729.737590Z_SC8_GSIR_RUBIN_SWF1_R1.fits"
            >>> fname.split('/')[-1]
            'EUC_SIM_LSST-R14S00-1492656-SCIENCE-z_07CD30AD0911-0042104_20201021T162729.737590Z_SC8_GSIR_RUBIN_SWF1_R1.fits'
            >>> name = fname.split('/')[-1]
            >>> name.split('-')[2]
            '1492656'
        """
        try:
            # extract only file name
            name = fname.split('/')[-1]
            # now extract id visit
            ret = name.split('-')[2]
        except:
            print("ERROR: fname empty ?")
            print(f"fname : '{fname}'")
            raise
        return ret
    
    def extract_all_visit_ids(self):
        """ return dictionary of list file with same id visit
        """        
        dict_visit_imag = {}
        for fimag in self.l_select_images:
            id_visit = self.extract_visit_ids_from_name(fimag)
            # print(id_visit)    
            if id_visit in dict_visit_imag:
                if len(dict_visit_imag[id_visit]) < self.limit_max_file:
                    dict_visit_imag[id_visit].append(fimag)
            else:
                dict_visit_imag[id_visit] = [fimag]
                print("add id visit " + id_visit)
        self.dict_visit_imag = dict_visit_imag
        self.nb_visit = len(dict_visit_imag)
    
    def create_link_butler(self):
        for idx , id_visit in enumerate(self.dict_visit_imag):
            lk_butler = os.path.join(self.params.root_butlers, self.dir_butler.replace('x', str(idx)))
            butler = os.path.join(self.params.root_butlers, self.dir_butler.replace('x', id_visit))
            try:          
                os.symlink(butler, lk_butler)
            except FileExistsError:
                pass            
    
    def search_max_image_by_visit(self):
        self.max_file = 0
        for idv, limag in self.dict_visit_imag.items():
            nb_image = len(limag)
            if nb_image > self.max_file:
                self.max_file = nb_image
        
    def create_file_image_by_visit(self):
        """ in other word it's the load balancing for array job
        """        
        lb_dir = os.path.join(self.params.root_butlers, "lb_image")
        try:
            os.mkdir(lb_dir)
        except FileExistsError:
            pass        
        # here idx is the rank in array job
        for idx , key in enumerate(self.dict_visit_imag):
            visit_file = os.path.join(lb_dir, self.file_lb.replace('xx', str(idx)))
            f_visit = open(visit_file, 'w')
            for fv in self.dict_visit_imag[key]:
                f_visit.write(fv + "\n")                
            f_visit.close()

    def select_raw_image(self):
        '''
        create list image
        '''
        lines = [line.strip('\n ') for line in open(self.file_name_image, 'r') if line.strip('\n ') != ""]
        if self.pattern:
            lines_visit = []
            for pattern in  self.pattern:                
                lines_visit += [line for line in lines if pattern in line]
        else:
            lines_visit = lines            
        self.l_select_images = lines_visit
                

    def get_dir_star_cat(self):
        return os.path.join(self.params.root_butlers, 'stars_cat')
        
    def select_and_loadbalancing(self):
        """ simple load balancing : image with same id visit are processing in same job 
        (because butler has no concurrent access in CCIN2P3)
        """
        # use filter file
        if not os.path.exists(self.params.root_butlers):
            os.system("mkdir " + self.params.root_butlers)
            os.system("mkdir " + self.params.root_butlers + '/log_job_lsst')
            os.system("mkdir " + self.params.root_butlers + '/log_job_euclid')       
            os.system("mkdir " + self.get_dir_star_cat())
            self.select_raw_image()
            self.extract_all_visit_ids()
            self.create_file_image_by_visit()     
            self.search_max_image_by_visit()
        else:
            # to initialize atttribut
            print("Skip butler creatio")
            self.select_raw_image()
            self.extract_all_visit_ids()
            self.create_file_image_by_visit()
            self.search_max_image_by_visit()
    
    def _line_file_to_list(self, n_file):         
        return [line.strip('\n ') for line in open(n_file, 'r') if line.strip('\n ') != ""]    
        
    def set_loadbalancing(self): 
        print("Read file ", self.file_name_image)       
        self.l_select_images = self._line_file_to_list(self.file_name_image)
    
    def retrieve_cat_name(self, id_visit):
        '''return path/name of stars catalog file or None associated to id_visit
        
        :param id_visit:
        '''
        id_visit_num = int(id_visit)
        command = f"grep '\-{id_visit_num}_' " + self.params.file_list_star_cat
        print(command)
        args = shlex.split(command)
        out_proc = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        output = out_proc.stdout.decode("utf-8").strip('\n ').split('\n')        
        print(output)
        if out_proc.returncode != 0:
            return None
        return output[0]
    
        
    def check_butler(self):
        '''
        * extract name of stars catalog associated to current id_visit
        *  create directory butler name : self.path_butler
        * launch script butler4euclid.py with this file if butler doesn't exist.
        '''
        id_visit = self.extract_visit_ids_from_name(self.l_select_images[0])        
        self.path_butler = os.path.join(self.params.root_butlers, self.dir_butler.replace('x', id_visit))
        if not self.stack.is_ref_cats(self.path_butler):
            file_cat = self.retrieve_cat_name(id_visit)
            if file_cat is None:
                return False            
            path_script = os.path.join(self.params.root_script, "at_ccin2p3")
            if False:
                command = f"{path_script}/butler4euclid.py -b {self.path_butler} -c {file_cat}"
                print(command)
                args = shlex.split(command)
                print(args)
                out_proc = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                output = out_proc.stdout        
                print(output.decode("utf-8"))
                if out_proc.returncode != 0:
                    return False
            o_but = butler.CreateButlerFromEuclidSimu(self.stack)
            o_but.init_file(self.path_butler, file_cat)
            o_but.create_insert_ref_cat()
            return o_but.status == 0            
        else:
            print('Butler already exists, skip creating, ingest ref cat ...')           
            return True            
    
                        
    def lsst_processing(self):
        """ main loop : 1)read load balancing, 2)create butler, 3) ingest image 4)process 5)extract FITS
        """
        # 1)
        self.set_loadbalancing()
        nb_image = len(self.l_select_images)
        print('\n-----------------------------------')
        print('{0} image(s) selected for processing:'.format(nb_image))
        print('-----------------------------------')
        for mfile in self.l_select_images:
            print(mfile)        
        cpt_error = 0
        print('')
        # 2) check butler
        if not self.check_butler():
            print("Can't create butler !!!")
            return 100        
        if self.params.del_fits:
            # remove all fits 
            mypath = self.path_butler + "/stage1"
            files_to_remove = glob.glob(os.path.join(mypath, '*.fits*'))
            for oldfile in files_to_remove:
                os.remove(oldfile)        
        
        for idx, mfile in enumerate(self.l_select_images):            
            print("============================= {0}/{1}".format(idx + 1, nb_image))
            self.abs_cur_image = mfile
            self.name_cur_image = os.path.basename(mfile)
            # 3) 4)            
            ret = self.ingest_process_image()
            cpt_error += ret
            if ret == 0:
                self.extract_fits_euclid_from_butler()
            # 5) TO DO
        print("\n\nTotal error LSST: ", cpt_error)
        return cpt_error

    def eas_ingestion_only(self):
        # 1) read load balancing
        lb_dir = os.path.join(self.params.root_butlers, "lb_eas")
        f_lb = os.path.join(lb_dir, f"lb_eas_{self.rank}.txt")
        l_sef = self._line_file_to_list(f_lb)
        
        #2) ingestion SEF one by one
        print(f"EAS ingestion of {len(l_sef)} SEF")
        path_ingest = self.params.eas_ingest_script
        user = self.params.user_eas
        path_pwd = self.params.file_pwd_eas
        EAS_pjt = self.params.eas_project
        sdc = self.params.sdc_eas
        nb_nok = 0
        for sef in l_sef:
            visit = extract_visit_from_xml_name(sef)
            path_butler = os.path.join(self.params.root_butlers,f"butler_{visit}/stage1")
            if not os.path.exists(path_butler):
                # patch bug : vistid  int or string
                path_butler = os.path.join(self.params.root_butlers,f"butler_0{visit}/stage1")
            os.chdir(path_butler)
            command = f"{path_ingest} store {sef}.xml --username={user} "
            command += f"--password={path_pwd} --environment=test --SDC={sdc} --project={EAS_pjt}"            
            command += " --useoldfiles"
            print("Command : ", command)
            args = shlex.split(command)
            out_proc = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            output = out_proc.stdout.decode("utf-8")
            print(output)
            if 'Error' in output:
                nb_nok += 1
        print("\n\nTotal error EAS: ", nb_nok)
        return nb_nok

    def euclid_processing(self):
        """ main loop : 1)read load balancing 2) create xml 3) ingest data/XML in EAS
        """
        cpt_error = 0
        # 1) read load balancing
        self.set_loadbalancing()
        self.check_butler()
        # 2) create xml
        if not self.params.skip_xml:            
            #  build.x86_64-co7-gcc48-o2g/run InstallArea/x86_64-co7-gcc48-o2g/scripts/script_createXML -p ~/Work/Projects/testfits
            os.chdir(f"{self.params.root_script}/EXT_PF1_dm4lsst")     
            path_run = "build.x86_64-conda_cos6-gcc73-o2g/run"
            path_script = f"{self.params.root_script}/EXT_PF1_dm4lsst/InstallArea/x86_64-conda_cos6-gcc73-o2g/scripts"
            command = f"{path_run} {path_script}/script_createXML -p {self.path_butler}/stage1 -t {self.params.pipeline_run}"
            if self.params.data_set_release != "":
                command += f" -r {self.params.data_set_release}"
            if self.params.filter_pos != "":
                command += f" -f {self.params.filter_pos}"
            print("Command : ", command)
            args = shlex.split(command)
            out_proc = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            output = out_proc.stdout
            print(output.decode("utf-8"))
            if out_proc.returncode != 0:
                print("ERROR : " + command)
                return 1
        # 3) ingest data/XML in EAS
        # ./ingest_client_sc456.py store $PWD/testfits  --username=jcolley --password=/home/user/Work/Projects/mypwd.txt --environment=current --SDC=SDC-FR --project=TEST
        if self.params.skip_ingest:
            print("Skip EAS ingestion")
            return 0
        path_ingest = self.params.eas_ingest_script    
        user = self.params.user_eas
        path_pwd = self.params.file_pwd_eas
        EAS_pjt = self.params.eas_project
        sdc = self.params.sdc_eas
        command = f"{path_ingest} store {self.path_butler}/stage1 --username={user} "
        command += f"--password={path_pwd} --environment=test --SDC={sdc} --project={EAS_pjt}"
        if self.params.use_old_files:
            command += " --useoldfiles"
        print("Command : ", command)
        args = shlex.split(command)
        out_proc = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        output = out_proc.stdout
        print(output.decode("utf-8"))
        if out_proc.returncode != 0:
            print("ERROR : " + command)
            return 1.
        print("\n\nTotal error Euclid: ", cpt_error)
        return cpt_error

    def ingest_process_image(self):
        """ change name image, insert in butler and process it
        """
        #  1/3 :change name and create a link        
        hdul = fits.open(self.abs_cur_image)
        name_lsst = hdul[0].header['OUTFILE']
        hdul.close()
        print('Name euclid :', self.name_cur_image)
        print('Name LSST   :', name_lsst)        
        dlsst = analyse_name_ccd_file(name_lsst)
        self.dlsst = dlsst
        print('path butler: ', self.path_butler)
        os.chdir(self.path_butler)
        lk_file = "lk_data/" + name_lsst
        if not os.path.exists(lk_file):            
            print('filter : change LSST_x to x')
            os.system('cp ' + self.abs_cur_image + ' ' + lk_file)
            hdul = fits.open(lk_file, 'update')
            filter_kw = hdul[0].header['FILTER']
            if filter_kw[:5] == "LSST_":
                hdul[0].header['FILTER'] = filter_kw[5].lower()
            hdul.flush()
            hdul.close()   
#             else:      
#                 os.symlink(pathfile, lk_file)
            #  2/3 : insert in butler
            command = self.stack.ingest_image_cmd(lk_file)
            print("Command : ", command)
            args = shlex.split(command)
            out_proc = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            output = out_proc.stdout
            print(output.decode("utf-8"))
            if out_proc.returncode != 0:
                print("ERROR : " + command)
                return 1.
            # Remove  temp raw image
            # os.system('rm '+lk_file)
        else:
            print("Skip insert image, already in butler !!!")
        # 3/3 : process image
        if self.skip_proc:
            print('Skip processing processEimage')
            return 0            
        command = self.stack.process_image_cmd(dlsst['visit'], dlsst['raft'], dlsst['sensor'])        
        print("Command : ", command)
        args = shlex.split(command)
        # out_proc = subprocess.run(args, env=self.lsst_env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)        
        if False:
            out_proc = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            output = out_proc.stdout
            print(output.decode("utf-8"))
            if out_proc.returncode != 0:
                print("ERROR : " + command)
                return 1
        with subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT) as out_proc:
#         out_proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
#         if True:
            max_mem = int(3.1 * 1024 ** 3)
            timer_s = 0.1 
            try:
                info_proc = psutil.Process(out_proc.pid)
                proc_exist = True
            except:
                # process already finished
                proc_exist = False
                
            while proc_exist:
                try:
                    # check memory
                    out_proc.poll()
                    # info_proc = psutil.Process(out_proc.pid)                    
                    rss = info_proc.memory_info().rss
                    # print(rss)
                    if rss < max_mem:
                        time.sleep(timer_s)
                    else:                        
                        rss_mb = int(rss / 1024 ** 2)
                        print(f"processEimage uses too many memory {rss_mb}MB")
                        print("kill exces mem: " + command)
                        # out_proc.kill()
                        out_proc.terminate()
                        break
                except:
                    # process dead
                    # print("finished")
                    break
            out_proc.wait()
            try:
                output = out_proc.stdout.read()
                print(output.decode("utf-8"))
            except:
                print("=====================PB output")                
            ret_proc = out_proc.returncode
            # print('ret proc:', ret_proc)
        if ret_proc != 0:
            print("ERROR : " + command)
            return 1
        return 0
    
    def extract_fits_euclid_from_butler(self):
        # write le fits detrended framepath_butler
        # python lsst2euc.py $BUTLER --visit 100176 --raft 0,1 --sensor 0,0 --force
        command = f"python {self.params.root_script}/convert_tools/lsst2euc.py {self.path_butler} --visit {self.dlsst['visit']}"        
        command += f" --raft {self.dlsst['raft']} --sensor {self.dlsst['sensor']} -o {self.path_butler}/stage1 --force -t {self.params.pipeline_run}"
        print("Command : ", command)
        args = shlex.split(command)        
        out_proc = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        output = out_proc.stdout
        print(output.decode("utf-8"))
        if out_proc.returncode != 0:
            print("ERROR : " + command)
            return 1        
        return 0
    
    def create_script_sge_for_lsst(self):
        script_sge = '''#!/bin/sh
#$ -l os=cl7 
#$ -l sps=1
#$ -P P_euclid_ext
#$ -N ext-stage1
'''
        # estimation of time processing for scheduler    
        cpu_time_image = 80
        if self.skip_proc:
            cpu_time_image += 100
        factor_safe = (3.2 * (self.max_file - 1) + (189 - self.max_file) * 3) / (189.0 - 1)
        print(f"factor_safe {factor_safe} for {self.max_file} image by proc.")
        butler_creation = 200
        extract_fits = 100
        tot_time = int(self.max_file * (cpu_time_image + extract_fits) * factor_safe + butler_creation)        
        script_sge += "#$ -l h_cpu=" + str(tot_time)
        script_sge += "\n#$ -t 1-" + str(self.nb_visit)
        script_sge += "\n#$ -wd " + self.params.root_butlers + "/log_job_lsst"
        script_sge += '''
#$ -l s_rss=3G
#$ -m be

# source LSST
source /cvmfs/sw.lsst.eu/linux-x86_64/lsst_distrib/w_2018_31/loadLSST.bash
setup lsst_distrib

# start main script loop on raw image lsst with some number visit
export OPENBLAS_NUM_THREADS=1
'''
        script_sge += f"cd {self.params.root_script}"
        script_sge += f"\n. init_stage1.bash"
        script_sge += "\nexport PATH_FP=" + self.params.nfile
        script_sge += '''

# lancement du job lsst
'''
        script_sge += "python3 at_ccin2p3/script_lsst_step.py $PATH_FP"        
        file_job = self.params.root_butlers + "/sge_arrayjob_lsst.sh"
        with open(file_job, 'w') as f_job:
            f_job.write(script_sge)
        print('Create ' + file_job + " , check it and do qsub " + file_job)

    def create_script_sge_for_eas(self, nb_sef_by_job):
        script_sge = '''#!/bin/sh
#$ -l os=cl7 
#$ -l sps=1
#$ -P P_euclid_ext
#$ -N eas-lsst
'''
        # 3600 mean by image
        # factor security : 2
        cpu_time_image = 60
        factor_safe = 1.2
        tot_time = int(nb_sef_by_job * cpu_time_image * factor_safe)        
        script_sge += "#$ -l h_cpu=" + str(tot_time)
        script_sge += "\n#$ -t 1-" + str(self.params.nb_proc_ingest)
        script_sge += "\n#$ -wd " + self.params.root_butlers + "/log_job_euclid"
        script_sge += '''
#$ -l s_rss=3G
#$ -m a

# source EDEN
source /cvmfs/euclid-dev.in2p3.fr/CentOS7/EDEN-2.1/bin/activate

'''
        script_sge += "export PATH_SCRIPT=" + os.path.dirname(__file__)             
        script_sge += "\nexport PATH=$PATH:$PATH_SCRIPT"
        script_sge += "\nexport PATH_FP=" + self.params.nfile
        script_sge += '''

# lancement du job euclid
'''
        script_sge += "python3 $PATH_SCRIPT/script_ingest_eas.py $PATH_FP"
        file_job = self.params.root_butlers + "/sge_arrayjob_eas.sh"
        with open(file_job, 'w') as f_job:
            f_job.write(script_sge)
        print('Create ' + file_job + " , check it and do qsub " + file_job)

    def create_script_sge_for_euclid(self):
        script_sge = '''#!/bin/sh
#$ -l os=cl7 
#$ -l sps=1
#$ -P P_euclid_ext
#$ -N ext-stage1
'''
        # 3600 mean by image
        # factor security : 2
        cpu_time_image = 60
        factor_safe = 1.2
        tot_time = int(self.max_file * cpu_time_image * factor_safe)        
        script_sge += "#$ -l h_cpu=" + str(tot_time)
        script_sge += "\n#$ -t 1-" + str(self.nb_visit)
        script_sge += "\n#$ -wd " + self.params.root_butlers + "/log_job_euclid"
        script_sge += '''
#$ -l s_rss=3G
#$ -m be

# source EDEN
source /cvmfs/euclid-dev.in2p3.fr/CentOS7/EDEN-2.1/bin/activate

'''
        script_sge += "export PATH_SCRIPT=" + self.params.root_script           
        script_sge += "\nexport PATH=$PATH:$PATH_SCRIPT"
        script_sge += "\nexport PATH_FP=" + self.params.nfile
        script_sge += '''

# lancement du job euclid
'''
        script_sge += "python3 $PATH_SCRIPT/at_ccin2p3/script_euclid_step.py $PATH_FP"
        file_job = self.params.root_butlers + "/sge_arrayjob_euclid.sh"
        with open(file_job, 'w') as f_job:
            f_job.write(script_sge)
        print('Create ' + file_job + " , check it and do qsub " + file_job)
        
    def create_script_sge_for_stage1_1step(self):
        script_sge = '''#!/bin/sh
#$ -l os=cl7 
#$ -l sps=1
#$ -P P_euclid_ext
#$ -N lsst-s1
'''
        # 3600 mean by image
        # factor security : 2
        # estimation of time processing for scheduler    
        cpu_time_image = 80
        if self.skip_proc:
            cpu_time_image += 100
        factor_safe = (3.2 * (self.max_file - 1) + (189 - self.max_file) * 3) / (189.0 - 1)
        print(f"factor_safe {factor_safe} for {self.max_file} image by proc.")
        butler_creation = 200
        extract_fits = 100
        tot_time_lsst = int(self.max_file * (cpu_time_image + extract_fits) * factor_safe + butler_creation)        
        
        cpu_time_image = 60
        factor_safe = 1.2
        tot_time_euclid = int(self.max_file * cpu_time_image * factor_safe)        
        script_sge += "#$ -l h_cpu=" + str(tot_time_lsst + tot_time_euclid)
        script_sge += "\n#$ -t 1-" + str(self.nb_visit)
        script_sge += "\n#$ -wd " + self.params.root_butlers + "/log_job_euclid"
        script_sge += '''
#$ -l h_rss=3500M
#$ -m a

# source EDEN
source /cvmfs/euclid-dev.in2p3.fr/CentOS7/EDEN-2.1/bin/activate

'''
        script_sge += "export PATH_SRC=" + self.params.root_script
        script_sge += "\nexport PATH_FP=" + self.params.nfile
        script_sge += f"\n\ncd {self.params.root_script}"
        script_sge += f"\n. init_stage1.bash"

        script_sge += '''

# step 1: processing with stack lsst under Singularity
'''
        script_sge += f"singularity exec --bind /sps,$HOME {self.params.singularity_lsst_stack} "
        script_sge += f"$PATH_SRC/singularity/script_step_lsst.bash $PATH_FP $PATH_SRC"       
        script_sge += '''
       
# step 2: xml creating and ingest in EAS
python3 $PATH_SRC/at_ccin2p3/script_euclid_step.py $PATH_FP
'''       
        file_job = self.params.root_butlers + "/sge_stage1_1step.sh"
        with open(file_job, 'w') as f_job:
            f_job.write(script_sge)
        print('Create ' + file_job + " , check it and do qsub " + file_job)
        

class ProcessStage1LsstEasIn(ProcessStage1LSST):
    
    def __init__(self):
        super().__init__()
    
    def select_raw_image(self):
        #
        eas_in = eas.EasToolsToDownload()
        eas_in.verbose = True
        eas_in.set_params(self.params)
        print('Retrieve raw images name from EAS.')
        l_raw = eas_in.query_list_raw_image(self.params.data_in_sps_or_eas) 
        if l_raw:
            nb_raw = len(l_raw)
            print(f'Write list raw images from EAS. {nb_raw} images to process for DataSetRelease {self.params.data_in_sps_or_eas}')
            write_file_from_list(l_raw, self.params.file_list_images)
        else:
            print("Error can't retreieve name file raw image with DataSetRelease {self.params.data_in_sps_or_eas}")
            os.system(self.params.file_name_images)
        super().select_raw_image()
        
    def retrieve_cat_name(self, id_visit):
        '''return path/name of stars catalog file or None associated to id_visit
        
        Q1: PipelineDefinitionId
            curl "https://eas-dps-cus.test.euclid.astro.rug.nl/COORS?class_name=DpdExtLssSingleVisitImage&
            Data.DataStorage.DataContainer.FileName.DataFileName=$rawFile&
            project=TEST&
            fields=Header.PipelineDefinitionId"
        Q2: Name file of stars catalog
            curl "https://eas-dps-cus.test.euclid.astro.rug.nl/COORS?class_name=DpdTrueUniverseOutput&
            Header.DataSetRelease=LSST_SC8_703_R1&Data.EuclidPointingId=154608&
            Header.PipelineDefinitionId=SIM_EXT-LSST&
            project=TEST&
            fields=Data.StarCatalogFitsFile.DataContainer.FileName.DataFileName"        
        Q3: Download stars catalog
            curl --netrc-file /pbs/home/c/colley/lpwd_EAS_curl.txt 
            "https://dss-mdb.euclid.astro.rug.nl/$CatFile" -o $CatFile

        :param id_visit:
        '''
        eas_in = eas.EasToolsToDownload()
        eas_in.verbose = True
        eas_in.set_params(self.params)
        # Q1
        raw_file = self.l_select_images[0]
        def_id = eas_in.query_definition_id(raw_file)
        if not def_id:
            print("Can't retrieve PipelineDefinitionId")
            return ""        
        print('PipelineDefinitionId : ',def_id)
        # Q2
        name_stars_cat = eas_in.query_name_stars_catalog(
                            self.params.data_in_sps_or_eas,
                            id_visit,
                            def_id)
        if not name_stars_cat:
            print("Can't retrieve PipelineDefinitionId")
            return ""        
        print('Stars cat: ', name_stars_cat)      
        # Q3
        path_cat = os.path.join(self.get_dir_star_cat(), name_stars_cat)
        eas_in.cmd_curl_download(name_stars_cat, self.get_dir_star_cat())
        return path_cat
        
                    
    def ingest_process_image(self):
        """
        add raw image download from EAS
        """
        eas_in = eas.EasToolsToDownload()
        eas_in.verbose = True
        eas_in.set_params(self.params)        
        dest_path = self.path_butler+'/lk_data'
        path_name =  os.path.join(dest_path, self.name_cur_image)
        self.abs_cur_image = path_name
        if not os.path.exists(path_name):
            eas_in.cmd_curl_download(self.name_cur_image, dest_path)
        return super().ingest_process_image()