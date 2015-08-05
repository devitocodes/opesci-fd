import os
pwd = os.getcwd()
data = [[]] * 9
for f in os.listdir(pwd):
	if f.endswith('.cpp'):
		cmd = 'icl -openmp ' + f
		os.system(cmd)
		exe = f.replace('.cpp', '.exe')
		txt = f.replace('.cpp', '.txt')
		cmd = exe + ' > ' + txt
		os.system(cmd)
		f2 = open(txt)
		
