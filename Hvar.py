import numpy as np



def get_r(x,y,z):
    return np.sqrt(x*x + y*y + z*z)

def get_theta(z,r):
    return np.arccos(z/r)

def get_phi(y,r,theta):
    return np.arcsin(y/(r*np.sin(theta)))


def to_polar(x,y,z):
    r = get_r(x,y,z)
    theta = get_theta(z,r)
    phi = get_phi(y,r,theta)
    return r,theta,phi

def to_cart(r,theta,phi):
    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)
    return x,y,z


# unit basis vectors for polar coordinates
# derivatives of positional vector w.r.t r,theta,and phi
def polar_basis_vectors(r,theta,phi):
    i_vec = np.array([1,0,0])
    j_vec = np.array([0,1,0])
    k_vec = np.array([0,0,1])
    e_r = np.sin(theta)*np.cos(phi)*i_vec +np.sin(theta)*np.sin(phi)*j_vec + np.cos(theta)*k_vec
    e_theta = np.cos(theta)*np.cos(phi)*i_vec + np.cos(theta)*np.sin(phi)*j_vec - np.sin(theta)*k_vec
    e_phi = -np.sin(phi)*i_vec + np.cos(phi)*j_vec # missing k basis vector
    return e_r, e_theta, e_phi


# Faki potential E = x^3y^2z
def potential_cart(x,y,z):
    return x**3*y**2*z

def force_cart(x,y,z):
    f_x = 3*x**2*y**2*z
    f_y = 2*x**3*y*z
    f_z = x**3*y**2
    return np.array([f_x,f_y,f_z])



def potential_polar(r,theta,phi):
    return r**(6)*np.sin(theta)**(5)*np.cos(theta)*np.cos(phi)**3*np.sin(phi)**2

def force_polar(r,theta,phi):
    e_r, e_theta, e_phi = polar_basis_vectors(r,theta,phi)
    sin_phi = np.sin(phi)
    cos_phi = np.cos(phi)
    sin_theta = np.sin(theta)
    cos_theta = np.cos(theta)
    # dr = 6*r**(5)*sin_theta**(5)*np.cos(theta)*np.cos(phi)**3*np.sin(phi)**2
    # dtheta = r**(6)*(5*np.sin(theta)**4*np.cos(theta)**2-np.sin(theta)**6)*np.cos(phi)**3*np.sin(phi)**2 
    # dphi = r**(6) * np.sin(theta)**(5) * np.cos(theta) * (3*np.cos(phi)**2 * np.sin(phi) * np.sin(phi)**2 - 2 * np.cos(phi)**(3) * np.sin(phi))
    dr = 6 * r**5 * sin_theta**5 * cos_theta * cos_phi**3 * sin_phi**2
    dtheta = r**6 * (5 * sin_theta**4 * cos_theta**2 - sin_theta**6) * cos_phi**3 * sin_phi**2
    dphi = r**6 * sin_theta**5 * cos_theta * (3 * cos_phi**2 * sin_phi**3 - 2 * cos_phi**3 * sin_phi)
    print("e_r",e_r,"e_theta", e_theta,"e_phi", e_phi)
    return e_r*dr+(1/r)*e_theta*dtheta+(1/(r*sin_theta))*e_phi*dphi
    
def numerical_grad(x,y,z):
    h =0.00001
    f_x = (potential_cart(x+h,y,z)-potential_cart(x-h,y,z))/(2*h)
    f_y = (potential_cart(x,y+h,z)-potential_cart(x,y-h,z))/(2*h)
    f_z = (potential_cart(x,y,z+h)-potential_cart(x,y,z-h))/(2*h)
    return np.array([f_x,f_y,f_z])


def numerical_grad_polar(x,theta,phi):
    h =0.00001
    e_r, e_theta, e_phi = polar_basis_vectors(r,theta,phi)
    f_r = (potential_polar(r+h,theta,phi)-potential_polar(r-h,theta,phi))/(2*h)
    f_theta = (potential_polar(r,theta+r*h,phi)-potential_polar(r,theta-r*h,phi))/(2*r*h)
    f_phi = (potential_polar(r,theta,phi+r*np.sin(theta)*h)-potential_polar(r,theta,phi-r*np.sin(theta)*h))/(2*r*np.sin(theta)*h)
    return f_r*e_r+f_theta*e_theta+f_phi*e_phi

def conversion_matrix(theta,phi):
    cm = np.array([[np.sin(theta)*np.cos(phi),np.cos(theta)*np.cos(phi),-np.sin(phi)],[np.sin(theta)*np.sin(phi),np.cos(theta)*np.sin(phi),np.cos(phi)],[np.cos(theta),-np.sin(theta),0]])
    return cm


if __name__== "__main__":
    x = 5.4
    y = 3.2
    z = 3.7    
    r,theta,phi = to_polar(x,y,z)
    print("potential cart: ",potential_cart(x,y,z))
    print("potential polar:",potential_polar(r,theta,phi))
    print("analytical forces cartesian:",force_cart(x,y,z))#/numerical_grad(x,y,z))
    print("analytical force polar:",force_polar(r,theta,phi))
    print("numerical force polar:",numerical_grad_polar(r,theta,phi))
    e_r, e_theta, e_phi = polar_basis_vectors(r,theta,phi)
    cm = conversion_matrix(theta,phi)
    test_vec = 5*e_r+2.3*e_theta+4.3*e_phi
    print(np.matmul(cm,test_vec))
    print(test_vec)
