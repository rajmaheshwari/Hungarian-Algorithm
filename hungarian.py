import numpy as np
import math
np.set_printoptions(precision=2)

"""Calculate distances between all grid points and given input points"""
"""Grid starts from (0,0) and ends at (a,a) and n=a*a (for square matrix)"""
def calculate_distance_matrix(coord):
	n = coord.shape[0]
	a = int(math.sqrt(n))
	
	grid_coord = np.zeros((n,2))
	distance = np.zeros((n,n))
	
	for i in range(a):
		for j in range(a):
			grid_coord[a*i+j][0] = i
			grid_coord[a*i+j][1] = j
		
	for i in range(n):
		x_grid = grid_coord[i][0]
		y_grid = grid_coord[i][1]
		for j in range(n):
			x_diff = (coord[j][0]-x_grid)
			y_diff = (coord[j][1]-y_grid)
			distance[i][j] = math.sqrt(x_diff*x_diff + y_diff*y_diff) 
			
	return distance

"""Subtracting row and column minimum from each of the rows and columns"""
def step1(cost_matrix):
	row_min = np.min(cost_matrix,1)
	b = (cost_matrix.transpose() - row_min.transpose()).transpose()
	#~ print row_min
	#~ print b
	#~ print 

	col_min = np.min(b,0)
	c = (b - col_min)
	#~ print col_min
	#~ print c
	#~ print
	return c


"""Finding minimum number of lines covering all zeros in the matrix"""
def step2(c):
	print 'Function step2() begins&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
	n = c.shape[0]
	f = (c==0)
	f = 1-f.astype(int)
	f_original = 1 - (c==0).astype(int)
	d = sum(f.transpose()==0)
	print c
	print f
	print d
	print 
	zero_assigned = []
	row_assigned = np.zeros(n).astype(int)
	
	flag = 0
	j=0
	print 'Finding minimum no of lines to cover all zeros'
	while (j<=n):
		if flag==0:
			j=j+1
		else:
			j=1	
		for i in range(n):
			flag = 0
			if d[i]==j and row_assigned[i]==0:
				flag = 1
				row_assigned[i] = 1
				z = np.where(f[i]==0)
				col_index = z[0][0]
				f[:,col_index] = 1
				f[i,:] = 1
				f[i,col_index] = 0
				zero_assigned.append((i,col_index))
				d = sum(f.transpose()==0)
				#~ print f
				print 'Rows assigned so far (1 corresponds to assigned row)'
				print row_assigned
				print
				break
	
	#~ print
	#~ print 'Assigned cell positions of zeros to rows in f matrix'
	#~ print zero_assigned			
	#~ print

	row_marked = np.zeros(n)
	col_marked = np.zeros(n)
		
	row_marked[np.where(row_assigned==0)[0]] = 1
	row_marked_indices = np.where(row_marked==1)[0]
	
	print '-----------------------------------------------------------------'
	
	marked_rows = list(row_marked_indices)
	row_marked_copy = np.copy(row_marked)
	col_marked_copy = np.copy(col_marked)	
		
	while True:
		for i in row_marked_indices:
			col_marked[np.where(f_original[i,:]==0)[0]] = 1
			col_marked_indices = np.where(col_marked==1)[0]
			#~ print str(col_marked_indices)+' are the columns marked'
			for j in col_marked_indices:
				for tup in zero_assigned:
					if tup[1]==j:
						marked_rows.append(tup[0])
						marked_rows = list(set(marked_rows))
				for k in marked_rows:
					row_marked[k] = 1	
				#~ print str(marked_rows)+' are the rows marked'
					
		if np.array_equal(row_marked_copy,row_marked) and np.array_equal(col_marked_copy,col_marked):
			break
		row_marked_copy = np.copy(row_marked)
		col_marked_copy = np.copy(col_marked)
	
		row_marked_indices = np.where(row_marked==1)[0]
	
	row_unmarked = 1 - row_marked
	print 'Row indices-'
	print row_unmarked
	print 'Column indices-'
	print col_marked_copy
	num_ass = sum(row_unmarked==1)+sum(col_marked_copy==1)
	print 'Function step2() ends&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
	return num_ass, row_unmarked, col_marked_copy

"""Finding the minimum value in matrix and subtracting and adding it"""
def step3(c,cols_marked,rows_marked):
	print 'Function step3() begins--------------------------------------'
	n = c.shape[0]
	temp = np.ones(c.shape)
	temp = np.multiply(temp,c)
	max_val = np.max(temp)
	
	col_indices = np.where(cols_marked==1)[0]
	row_indices = np.where(rows_marked==1)[0]
	
	for i in col_indices:
		temp[:,int(i)] = max_val
	for i in row_indices:
		temp[int(i),:] = max_val
	
	min_val = np.min(temp)
	
	for i in range(n):
		if i not in row_indices:
			c[int(i),:] = c[int(i),:] - min_val
	for i in col_indices:
		c[:,int(i)] = c[:,int(i)] + min_val		
	print 'Function step3() ends----------------------------------------'
	
	return c

"""Finding the final matching using the matrix obtained from applying above defined functions on the input"""
def step4(c):
	print 'Function step4() begins***********************************************'															
	n = c.shape[0]
	f = (c==0)
	f = 1-f.astype(int)
	f_original = 1 - (c==0).astype(int)
	d = sum(f.transpose()==0)
	
	zero_assigned = []
	row_assigned = np.zeros(n).astype(int)
	
	flag = 0
	j=0
	while (j<=n):
		if flag==0:
			j=j+1
		else:
			j=1	
		for i in range(n):
			flag = 0
			if d[i]==j and row_assigned[i]==0:
				flag = 1
				row_assigned[i] = 1
				z = np.where(f[i]==0)
				col_index = z[0][0]
				f[:,col_index] = 1
				f[i,:] = 1
				f[i,col_index] = 0
				zero_assigned.append((i,col_index))
				d = sum(f.transpose()==0)
				break
	
	print 'Function step4() ends***********************************************'
	return zero_assigned
	
"""Combining all the above functions and applying them until a solution is obtained"""
def hungarian(cost_matrix):
	count = 5
	c = step1(cost_matrix)
	while True:
		num_ass, final_rows_marked, final_cols_marked = step2(c)		
		if num_ass==c.shape[0]:
			break
		c = step3(c,final_cols_marked,final_rows_marked)
		count = count-1
		if count==0:
			print "This is the not good==============================================="
			break
			

	cells = step4(c)
	return cells

"""Maps input coordinates to the grid coordinates using the output of the method hungarian()"""
def map_points(cells,coord):
	a = int(math.sqrt(coord.shape[0]))
	print '----------------------------------------------'
	print
	print 'Input coordinates'
	print coord
	print
	for tup in cells:
		y = tup[0]%a
		x = (tup[0]-y)/a
		print 'Point '+str(coord[tup[1]])+' snaps to grid point '+str((x,y))+'.'
#~ ---------------------------------------------------------------------

n = 2
coord = n*np.random.rand(n*n,2)
cost_matrix = np.asarray(np.matrix([[ 0.85531522,  1.29252846],[ 1.93898172,  0.27637576],[ 0.55939688,  0.46164618],[ 0.51974249,  1.81025693],[ 0.50725915,  1.1675708 ],[ 0.84961689,  0.01215642], [ 1.27373429,  0.10067068],[ 1.62659856,  0.88620218],[ 1.35347636,  1.01616913]]))
cost_matrix = calculate_distance_matrix(coord)
cells = hungarian(cost_matrix)
map_points(cells,coord)
print cost_matrix
