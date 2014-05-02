nu = 1
tau = 1
h_x = 1
h_y = 1
N_X = 3
N_Y = 3
i_0 = 2
j_0 = 2
i_N_X = i_0 + N_X
j_N_Y = j_0 + N_Y
border = 0
omega = 1.65
maxiter = 20
class Cell:
    ''' ucc,ulc,urc,vcc,vcu,vcb'''
    def __init__(self,data=[]):
        self.data = data
        
def calc_DU_horizontal(left_point,point,right_point,second_right_point):
    '''u flux through !!!vertical!!! borders; point = (i,j)'''
    ucc = 0.5*(right_point.data[2] - point.data[2])
    delta_x = right_point.data[2] - point.data[2]
    delta_sqr_x = 0.5*(second_right_point.data[2]-right_point.data[2]-point.data[2]+left_point.data[2])
    if ucc*delta_x*delta_sqr_x >= 0:
        velocity_of_subst_transfer = MOR_U_horizontal(left_point,point,right_point,second_right_point)
    else:
        velocity_of_subst_transfer = MCR_U_horizontal(point,right_point)
    return ucc*velocity_of_subst_transfer

def MOR_U_horizontal(left_point,point,right_point,second_right_point):
    ucc = 0.5*(right_point.data[2] - point.data[2])
    c_point = tau*abs(ucc)/h_x
    if ucc >= 0:
        velocity_of_subst_transfer = 0.5*(3-c_point)*point.data[2]-0.5*(1-c_point)*left_point.data[2]
    else:
        velocity_of_subst_transfer = 0.5*(3-c_point)*right_point.data[2]-0.5*(1-c_point)*second_right_point.data[2]
    return velocity_of_subst_transfer

def MCR_U_horizontal(point,right_point):
    ucc = 0.5*(right_point.data[2] - point.data[2])
    velocity_of_subst_transfer = 0.5*(1-tau*ucc/h_x)*right_point.data[2]+0.5*(1+tau*ucc/h_x)*point.data[2]
    return velocity_of_subst_transfer


def calc_DU_vertical(bottom_point,point,upper_point,second_upper_point,right_point):
    '''u flux through !!!horizontal" borders; point = (i,j)'''
    vru = 0.5*(right_point.data[4]+point.data[4])
    delta_y = upper_point.data[2]-point.data[2]
    delta_sqr_y = second_upper_point.data[2] - upper_point.data[2] - point.data[2] + bottom_point.data[2]
    if vru*delta_y*delta_sqr_y >= 0:
        vost = MOR_U_vertical(bottom_point,point,upper_point,second_upper_point,right_point)
    else:
        vost = MCR_U_vertical(point,upper_point,right_point)
    return vost*vru

def MOR_U_vertical(bottom_point,point,upper_point,second_upper_point,right_point):
    vru = 0.5*(right_point.data[4]+point.data[4])
    point_c = abs(vru)*tau/h_y
    ''' WTF in da book'''
    if vru >= 0:
        velocity_of_subst_transfer = 0.5*(3-point_c)*point.data[2]-0.5*(1-point_c)*bottom_point.data[2]
    else:
        velocity_of_subst_transfer = 0.5*(3-point_c)*upper_point.data[2]-0.5*(1-point_c)*second_upper_point.data[2]
    return velocity_of_subst_transfer

def MCR_U_vertical(point,upper_point,right_point):
    vru = 0.5*(right_point.data[4]+point.data[4])
    velocity_of_subst_transfer = 0.5*(1-tau*vru/h_y)*upper_point.data[2]+0.5*(1+tau*vru/h_y)*point.data[2]
    return velocity_of_subst_transfer

def calc_DV_horizontal(left_point,point,right_point,second_right_point,upper_point):
    ''' written relatively to (i,j)'''
    uru = 0.5*(upper_point.data[2]-point.data[2])
    delta_x = right_point.data[4] - point.data[4]
    delta_sqr_x = second_right_point.data[4]-right_point.data[4]-point.data[4]+left_point.data[4]
    if uru*delta_x*delta_sqr_x >=0:
        velocity_of_subst_transfer = MOR_V_horizontal(left_point,point,right_point,second_right_point,upper_point)
    else:
        velocity_of_subst_transfer = MCR_V_horizontal(point,right_point,upper_point)
    return velocity_of_subst_transfer*uru
      
def MOR_V_horizontal(left_point,point,right_point,second_right_point,upper_point):
    uru = 0.5*(upper_point.data[2]-point.data[2])
    c_point = tau*abs(uru)/h_x
    if uru >= 0:
        velocity_of_subst_transfer = 0.5*(3-c_point)*point.data[4]-0.5*(1-c_point)*left_point.data[4]
    else:
        velocity_of_subst_transfer = 0.5*(3-c_point)*right_point.data[4]-0.5*(1-c_point)*second_right_point.data[4]
    return velocity_of_subst_transfer
        
def MCR_V_horizontal(point,right_point,upper_point):
    uru = 0.5*(upper_point.data[2]-point.data[2])
    velocity_of_subst_transfer = 0.5*(1-tau*uru/h_x)*right_point.data[4]+0.5*(1+tau*uru/h_x)*point.data[4]
    return velocity_of_subst_transfer

def calc_DV_vertical(bottom_point,point,upper_point,second_upper_point):
    ''' written relatively to (i,j)'''
    vcc = 0.5*(upper_point.data[4] + point.data[4])
    delta_y = upper_point.data[4] - point.data[4]
    delta_sqr_y = second_upper_point.data[4] - upper_point.data[4] - point.data[4] + bottom_point.data[4]
    if vcc*delta_y*delta_sqr_y >= 0:
        velocity_of_subst_transfer = MOR_V_vertical(bottom_point,point,upper_point,second_upper_point)
    else:
        velocity_of_subst_transfer = MCR_V_vertical(point,upper_point)
    return vcc*velocity_of_subst_transfer
    
def MOR_V_vertical(bottom_point,point,upper_point,second_upper_point):
    vcc = 0.5*(upper_point.data[4] - point.data[4])
    c_point = tau*abs(vcc)/h_y
    if vcc >= 0:
        velocity_of_subst_transfer = 0.5*(3-c_point)*point.data[4]-0.5*(1-c_point)*bottom_point.data[4]
    else:
        velocity_of_subst_transfer = 0.5*(3-c_point)*upper_point.data[4]-0.5*(1-c_point)*second_upper_point.data[4]
    return velocity_of_subst_transfer
    
def MCR_V_vertical(point,upper_point):
    vcc = 0.5*(upper_point.data[4] - point.data[4])
    velocity_of_subst_transfer = 0.5*(1-tau*vcc/h_y)*upper_point.data[4]+0.5*(1+tau*vcc/h_y)*point.data[4]
    return velocity_of_subst_transfer

def set_P_boundary(solution):
    '''bottom and upper'''
    for i in range(i_0,i_N_X):
        solution[i][j_0-1].data[5] = solution[i][j_0].data[6]
        solution[i][j_N_Y].data[6] = solution[i][j_N_Y-1].data[6]
    '''left and right'''
    for j in range(j_0,j_N_Y):
        solution[i_0-1][j].data[5] = solution[i_0][j].data[6]
        solution[i_N_X][j].data[6] = solution[i_N_X-1][j].data[6]
        
def calc_P(left_point,point,right_point,upper_point,bottom_point):
    point.data[5] = ((1-omega)*point.data[6]+omega*((1/h_x**2)*left_point.data[5]+(1/h_y**2)*bottom_point.data[5])/(2/h_x**2+2/h_y**2)
                     +((1/h_x**2)*right_point.data[6]+(1/h_y**2)*upper_point.data[6])/(2/h_x**2+2/h_y**2)
                     -omega/(2/h_x**2+2/h_y**2)*((point.data[1]-left_point.data[1])/h_x+(point.data[3]-bottom_point.data[3])/h_y)/tau)

def calc_new_velocity(point,right_point,upper_point):
    point.data[2] = point.data[1]-tau/h_x*(right_point.data[5]-point.data[5])
    point.data[4] = point.data[3]-tau/h_y*(upper_point.data[5]-point.data[5])
    
def update_P(solution):
    for i in range(i_0,i_N_X):
        for j in range(j_0,j_N_Y):
            solution[i][j].data[6] = solution[i][j].data[5]
            
def set_velocity_boundary(solution):
    '''v bottom wall'''
    v_bottom = 0
    for i in range(i_0,i_N_X):
        solution[i][j_0-1].data[4] = v_bottom
        for j in range(0,j_0-1):
            solution[i][j].data[4] = -solution[i][2*j_0-j-1].data[4]+v_bottom
    '''v upper wall'''
    v_upper = 0
    for i in range(i_0,i_N_X):
        solution[i][j_N_Y-1].data[4] = v_upper
        for j in range(0,j_0):
            solution[i][j_N_Y+j].data[4] = -solution[i][j_N_Y-2-j].data[4]
    
    '''right wall'''
    v_right = 0
    for j in range(j_0,j_N_Y):
        for i in range(0,i_0):
            solution[i_N_X+i][j].data[4] = -solution[i_N_X-1-i][j].data[4]+v_right
    u_right = 0
    for j in range(j_0,j_N_Y):
        solution[i_N_X-1][j].data[2] = u_right
        for i in range(0,i_0):
            solution[i_N_X+i][j].data[2] = -solution[i_N_X-2-i][j].data[2]+u_right
    
    '''left wall'''
    v_left = 0
    for j in range(j_0,j_N_Y):
        for i in range(0,i_0):
            solution[i][j].data[4] = -solution[2*i_0-i-1][j].data[4]+v_left
    u_left = 0
    for j in range(j_0,j_N_Y):
        solution[i_0-1][j].data[2] = u_left
        for i in range(0,i_0-1):
            solution[i][j].data[2] = -solution[2*i_0-1-i][j].data[2]+u_left
    '''u bottom wall'''
    u_bottom = 0
    for i in range(i_0,i_N_X):
        for j in range(0,j_0):
            solution[i][j].data[2] = -solution[i][2*j_0-j-1].data[2]+u_bottom
    
    '''u upper wall'''
    u_upper = 1
    for i in range(i_0,i_N_X):
        for j in range(0,j_0):
            solution[i][j_N_Y+j].data[2] = -solution[i][j_N_Y-1-j].data[2]+u_upper
            
def calc_laplas4V(left_point,left_upper_point,point,upper_point,right_point):
    laplas = (((right_point.data[4] - point.data[4])/h_x - (upper_point.data[2] - point.data[2])/h_y)
              - ((point.data[4] - left_point.data[4])/h_x - (left_upper_point.data[2] - left_point.data[2])/h_y))
    return laplas

def calc_laplas4U(bottom_point,point,upper_point,right_point,right_bottom_point):
    laplas = (((right_point.data[4] - point.data[4])/h_x - (upper_point.data[2] - point.data[2])/h_y)
              - ((right_bottom_point.data[4] - bottom_point.data[4])/h_x - (point.data[2] - bottom_point.data[2])))
    return laplas
        
        
            
            

if __name__ == '__main__':
    
    elem = 4
    
    solution1 = [[Cell([0,111,2222,333,4444,555,666]) for j in range(7)] for i in range(7)]
    for j in range(j_0,j_N_Y):
        for i in range(i_0,i_N_X):
            solution1[i][j].data[elem] = 3*(j-j_0)+(i-i_0)+1
            
    for j in range(j_0,j_N_Y):
        for i in range(i_0,i_N_X):
            solution1[i][j].data[elem-2] = (j-j_0)+3*(i-i_0)+1
            
    for j in reversed(range(0,j_N_Y+j_0)):
        for i in range(0,i_N_X+i_0):
            print("%3.0d" % solution1[i][j].data[elem],end=' ')
        print()
    print()
   
    
    '''Part1'''
    set_velocity_boundary(solution1)
    
    for j in reversed(range(0,j_N_Y+j_0)):
        for i in range(0,i_N_X+i_0):
            print("%3.0d" % solution1[i][j].data[elem],end=' ')
        print()
    print()
    
    for i in range(i_0+border,i_N_X-border):
        for j in range(j_0+border,j_N_Y-border):
            DU3 = calc_DU_horizontal(solution1[i-1][j],solution1[i][j],solution1[i+1][j],solution1[i+2][j])
            DU1 = calc_DU_horizontal(solution1[i-2][j],solution1[i-1][j],solution1[i][j],solution1[i+1][j])
            DU4 = calc_DU_vertical(solution1[i][j-1],solution1[i][j],solution1[i][j+1],solution1[i][j+2],solution1[i+1][j])
            DU2 = calc_DU_vertical(solution1[i][j-2],solution1[i][j-1],solution1[i][j],solution1[i][j+1],solution1[i+1][j-1])
            
            DV3 = calc_DV_horizontal(solution1[i-1][j],solution1[i][j],solution1[i+1][j],solution1[i+2][j],solution1[i][j+1])
            DV1 = calc_DV_horizontal(solution1[i-2][j],solution1[i-1][j],solution1[i][j],solution1[i+1][j],solution1[i-1][j+1])
            DV4 = calc_DV_vertical(solution1[i][j-1],solution1[i][j],solution1[i][j+1],solution1[i][j+2])
            DV2 = calc_DV_vertical(solution1[i][j-2],solution1[i][j-1],solution1[i][j],solution1[i][j+1])
            
            laplas4U = calc_laplas4U(solution1[i][j-1],solution1[i][j],solution1[i][j+1],solution1[i+1][j],solution1[i+1][j-1])
            solution1[i][j].data[1] = solution1[i][j].data[2] - tau/h_x*(DU3-DU1) - tau/h_y*(DU4-DU2) - nu*tau/h_y*(laplas4U)
            
            laplas4V = calc_laplas4V(solution1[i-1][j],solution1[i-1][j+1],solution1[i][j],solution1[i][j+1],solution1[i+1][j])
            solution1[i][j].data[3] = solution1[i][j].data[4] - tau/h_y*(DV3-DV1) - tau/h_x*(DV4-DV2) - nu*tau/h_x*(laplas4V)
            
    for j in reversed(range(0,j_N_Y+j_0)):
        for i in range(0,i_N_X+i_0):
            print("%3.0d" % solution1[i][j].data[3],end=' ')
        print()        
    """   
    '''Part 2'''
    for _ in range(maxiter):
        set_P_boundary(solution1) 
        for i in range(i_0,i_N_X):
            for j in range(j_0,j_N_Y):
                calc_P(solution1[i-1][j],solution1[i][j],solution1[i+1][j],solution1[i][j+1],solution1[i][j-1])
        update_P(solution1)
        print("Iteration %d done" %(_))
    '''Part 3'''
    
    for i in range(i_0,i_N_X):
        for j in range(j_0,j_N_Y):
            calc_new_velocity(solution1[i][j])
    """
