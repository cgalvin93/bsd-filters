#time python ~/desktop/run_scores.py

import os
'''
#for running IAM on remaining structures of fd set
os.chdir('3dk9')
os.system('~/desktop/Rosetta/main/source/bin/rosetta_scripts.default.macosclangrelease -l pdblist.txt -extra_res_fa fad.params -parser:protocol ~/desktop/bsd_filters/iam.xml -parser:view -out:file:score_only 3dk9-iam-score.sc')
os.chdir('..')
os.chdir('3dlc')
os.system('~/desktop/Rosetta/main/source/bin/rosetta_scripts.default.macosclangrelease -l pdblist.txt -extra_res_fa sam.params -parser:protocol ~/desktop/bsd_filters/iam.xml -parser:view -out:file:score_only 3dlc-iam-score.sc')
os.chdir('..')
os.chdir('3r2q')
os.system('~/desktop/Rosetta/main/source/bin/rosetta_scripts.default.macosclangrelease -l pdblist.txt -extra_res_fa gsh.params -parser:protocol ~/desktop/bsd_filters/iam.xml -parser:view -out:file:score_only 3r2q-iam-score.sc')
os.chdir('..')
'''
#for running all merttric scoring script on all the structures
#cd desktop/prjk/analysis/strc = STARTING DIRECTORY, contains:
#bre, cm, fd
#which each contain: 1f4p  1zk4  2ij2  2xbn	3dk9  3dlc  3r2q  3s6f
#of which 1f4p  1zk4   2xbn	3dk9  3dlc  3r2q  are working

os.chdir('bre')
os.chdir('1f4p')
os.system('~/desktop/Rosetta/main/source/bin/rosetta_scripts.default.macosclangrelease -l pdblist.txt -extra_res_fa fmn.params -parser:protocol ~/desktop/bsd_filters/all_metrics.xml -parser:view -out:file:score_only bre-1f4p-all-score.sc')
os.chdir('..')
os.chdir('1zk4')
os.system('~/desktop/Rosetta/main/source/bin/rosetta_scripts.default.macosclangrelease -l pdblist.txt -extra_res_fa nap.params -parser:protocol ~/desktop/bsd_filters/all_metrics.xml -parser:view -out:file:score_only bre-1zk4-all-score.sc')
os.chdir('..')
os.chdir('2xbn')
os.system('~/desktop/Rosetta/main/source/bin/rosetta_scripts.default.macosclangrelease -l pdblist.txt -extra_res_fa pmp.params -parser:protocol ~/desktop/bsd_filters/all_metrics.xml -parser:view -out:file:score_only bre-2xbn-all-score.sc')
os.chdir('..')
os.chdir('3dk9')
os.system('~/desktop/Rosetta/main/source/bin/rosetta_scripts.default.macosclangrelease -l pdblist.txt -extra_res_fa fad.params -parser:protocol ~/desktop/bsd_filters/all_metrics.xml -parser:view -out:file:score_only bre-3dk9-all-score.sc')
os.chdir('..')
os.chdir('3dlc')
os.system('~/desktop/Rosetta/main/source/bin/rosetta_scripts.default.macosclangrelease -l pdblist.txt -extra_res_fa sam.params -parser:protocol ~/desktop/bsd_filters/all_metrics.xml -parser:view -out:file:score_only bre-3dlc-all-score.sc')
os.chdir('..')
os.chdir('3r2q')
os.system('~/desktop/Rosetta/main/source/bin/rosetta_scripts.default.macosclangrelease -l pdblist.txt -extra_res_fa gsh.params -parser:protocol ~/desktop/bsd_filters/all_metrics.xml -parser:view -out:file:score_only bre-3r2q-all-score.sc')
os.chdir('..')
os.chdir('..')

os.chdir('cm')
os.chdir('1f4p')
os.system('~/desktop/Rosetta/main/source/bin/rosetta_scripts.default.macosclangrelease -l pdblist.txt -extra_res_fa fmn.params -parser:protocol ~/desktop/bsd_filters/all_metrics.xml -parser:view -out:file:score_only cm-1f4p-all-score.sc')
os.chdir('..')
os.chdir('1zk4')
os.system('~/desktop/Rosetta/main/source/bin/rosetta_scripts.default.macosclangrelease -l pdblist.txt -extra_res_fa nap.params -parser:protocol ~/desktop/bsd_filters/all_metrics.xml -parser:view -out:file:score_only cm-1zk4-all-score.sc')
os.chdir('..')
os.chdir('2xbn')
os.system('~/desktop/Rosetta/main/source/bin/rosetta_scripts.default.macosclangrelease -l pdblist.txt -extra_res_fa pmp.params -parser:protocol ~/desktop/bsd_filters/all_metrics.xml -parser:view -out:file:score_only cm-2xbn-all-score.sc')
os.chdir('..')
os.chdir('3dk9')
os.system('~/desktop/Rosetta/main/source/bin/rosetta_scripts.default.macosclangrelease -l pdblist.txt -extra_res_fa fad.params -parser:protocol ~/desktop/bsd_filters/all_metrics.xml -parser:view -out:file:score_only cm-3dk9-all-score.sc')
os.chdir('..')
os.chdir('3dlc')
os.system('~/desktop/Rosetta/main/source/bin/rosetta_scripts.default.macosclangrelease -l pdblist.txt -extra_res_fa sam.params -parser:protocol ~/desktop/bsd_filters/all_metrics.xml -parser:view -out:file:score_only cm-3dlc-all-score.sc')
os.chdir('..')
os.chdir('3r2q')
os.system('~/desktop/Rosetta/main/source/bin/rosetta_scripts.default.macosclangrelease -l pdblist.txt -extra_res_fa gsh.params -parser:protocol ~/desktop/bsd_filters/all_metrics.xml -parser:view -out:file:score_only cm-3r2q-all-score.sc')
os.chdir('..')
os.chdir('..')



os.chdir('fd')
os.chdir('1f4p')
os.system('~/desktop/Rosetta/main/source/bin/rosetta_scripts.default.macosclangrelease -l pdblist.txt -extra_res_fa fmn.params -parser:protocol ~/desktop/bsd_filters/all_metrics.xml -parser:view -out:file:score_only fd-1f4p-all-score.sc')
os.chdir('..')
os.chdir('1zk4')
os.system('~/desktop/Rosetta/main/source/bin/rosetta_scripts.default.macosclangrelease -l pdblist.txt -extra_res_fa nap.params -parser:protocol ~/desktop/bsd_filters/all_metrics.xml -parser:view -out:file:score_only fd-1zk4-all-score.sc')
os.chdir('..')
os.chdir('2xbn')
os.system('~/desktop/Rosetta/main/source/bin/rosetta_scripts.default.macosclangrelease -l pdblist.txt -extra_res_fa pmp.params -parser:protocol ~/desktop/bsd_filters/all_metrics.xml -parser:view -out:file:score_only fd-2xbn-all-score.sc')
os.chdir('..')
os.chdir('3dk9')
os.system('~/desktop/Rosetta/main/source/bin/rosetta_scripts.default.macosclangrelease -l pdblist.txt -extra_res_fa fad.params -parser:protocol ~/desktop/bsd_filters/all_metrics.xml -parser:view -out:file:score_only fd-3dk9-all-score.sc')
os.chdir('..')
os.chdir('3dlc')
os.system('~/desktop/Rosetta/main/source/bin/rosetta_scripts.default.macosclangrelease -l pdblist.txt -extra_res_fa sam.params -parser:protocol ~/desktop/bsd_filters/all_metrics.xml -parser:view -out:file:score_only fd-3dlc-all-score.sc')
os.chdir('..')
os.chdir('3r2q')
os.system('~/desktop/Rosetta/main/source/bin/rosetta_scripts.default.macosclangrelease -l pdblist.txt -extra_res_fa gsh.params -parser:protocol ~/desktop/bsd_filters/all_metrics.xml -parser:view -out:file:score_only fd-3r2q-all-score.sc')
