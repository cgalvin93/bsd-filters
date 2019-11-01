#time python ~/desktop/run_scores.py

import os

os.chdir('3dk9')
os.system('~/desktop/Rosetta/main/source/bin/rosetta_scripts.default.macosclangrelease -l pdblist.txt -extra_res_fa fad.params -parser:protocol ~/desktop/bsd_filters/iam.xml -parser:view -out:file:score_only 3dk9-iam-score.sc')
os.chdir('..')
os.chdir('3dlc')
os.system('~/desktop/Rosetta/main/source/bin/rosetta_scripts.default.macosclangrelease -l pdblist.txt -extra_res_fa sam.params -parser:protocol ~/desktop/bsd_filters/iam.xml -parser:view -out:file:score_only 3dlc-iam-score.sc')
os.chdir('..')
os.chdir('3r2q')
os.system('~/desktop/Rosetta/main/source/bin/rosetta_scripts.default.macosclangrelease -l pdblist.txt -extra_res_fa gsh.params -parser:protocol ~/desktop/bsd_filters/iam.xml -parser:view -out:file:score_only 3r2q-iam-score.sc')
os.chdir('..')
