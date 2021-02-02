/* A fragmentation module for the N-body integrator REBOUND by Anna Childs.
    Copyright (C) 2021 Anna Childs

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

const double min_frag_mass = 1.4e-8;
const double mass_limit = 2.8e-4;  //maximum mass of bodies that are allowed to fragment


#define MIN(a, b) ((a) > (b) ? (b) : (a))    ///< Returns the minimum of a and b
#define MAX(a, b) ((a) > (b) ? (a) : (b))    // Returns the maximum of a and b



struct collision_params
{
    int target;
    int projectile;
    double dx;
    double dy;
    double dz;
    double b;
    double Vix;
    double Viy;
    double Viz;
    double Vi;
    double l;
    double rho1;
    double rho;
    double cstar;
    double mubar;
    double mu;
    double QR;
    double QpRD;
    double V_esc;
    double separation_distance;
    double Mlr;
    double Mslr;
    double Q;
    double Mlr_dag;
    double Q_star;
    double vrel;
    double xrel;
}; 


void make_vector(double x1, double y1, double z1, double x2, double y2, double z2, double *x, double*y, double*z){   //Galilean transform
    *x = x2-x1;
    *y = y2-y1;
    *z = z2-z1;
}

double get_dot(double x1, double y1, double z1, double x2, double y2, double z2){ 
    return (x1*x2)+(y1*y2)+(z1*z2);
        }  //return dot product of two vectors

double get_mag(double x, double y, double z){
    return sqrt(pow(x,2)+pow(y,2)+pow(z,2));
        }   //return magnitude of vector

double hill_radius(struct reb_simulation* const r, double x, double y, double z, double mass){
    struct reb_particle* particles = r->particles;
    return get_mag(x,y,z)*pow(mass/(3*particles[0].m), 1./3.);
}




int tot_no_frags = 0;
void add_fragments(struct reb_simulation* const r, struct reb_collision c, struct collision_params *params){
    struct reb_particle* target = &(r->particles[params->target]);
    struct reb_particle* projectile = &(r->particles[params->projectile]);
    struct reb_particle com = reb_get_com_of_pair(*target, *projectile);

    double initial_mass = target -> m + projectile -> m;

    double remaining_mass = initial_mass - params->Mlr;
    
    double big_frags = 0;
    if (params->Mslr > 0){
        remaining_mass = remaining_mass -  params->Mslr;
        big_frags = 1;
    }

    int no_frags = remaining_mass/min_frag_mass;  //fragments are broken up into equal sizes
    double frag_mass = remaining_mass/no_frags;

    int new_bodies = no_frags + big_frags;

    double rtot = target -> r + projectile -> r;

    char hash[10];
    double mxsum[3] = {0,0,0};
    double mvsum[3] = {0,0,0};
    //target gets mass of Mlr and is assigned COM position and velocity;

    target -> m = params->Mlr;
    target->x = com.x;
    target->y = com.y;
    target->z = com.z;

    target->vx = com.vx;
    target->vy = com.vy;
    target->vz = com.vz;
    target->lastcollision = r->t;

    mxsum[0] = mxsum[0] + target->m*target->x;
    mxsum[1] = mxsum[1] + target->m*target->y;  
    mxsum[2] = mxsum[2] + target->m*target->z;

    mvsum[0] = mvsum[0] + target->m*target->vx;
    mvsum[1] = mvsum[1] + target->m*target->vy; 
    mvsum[2] = mvsum[2] + target->m*target->vz;         
    
    double theta_inc = (2*M_PI)/new_bodies;
    

    double unit_vix, unit_viy, unit_viz, zx, zy, zz, z, ox, oy, oz, o;

    unit_vix = params->Vix/params->vrel;  //unit vector parallel to target velocity
    unit_viy = params->Viy/params->vrel;
    unit_viz = params->Viz/params->vrel;

    zx = (params->Viy*params->dz - params->Viz*params->dy);                     // vector normal to the collision plane; vrel cross xrel
    zy = (params->Viz*params->dx - params->Vix*params->dz);
    zz = (params->Vix*params->dy - params->Viy*params->dx);

    z = get_mag(zx, zy, zz);

    zx = zx/z;          //unit vector
    zy = zy/z;
    zz = zz/z;


    ox = (zy*params->Viz - zz*params->Viy);                   // vector normal to target velocity in collision plane; z cross vrel
    oy = (zz*params->Vix - zx*params->Viz);
    oz = (zx*params->Viy - zy*params->Vix);

    o = get_mag(ox, oy, oz);

    ox = ox/o;      //unit vector
    oy = oy/o;
    oz = oz/o;


    double fragment_velocity = sqrt(1.1*pow(params->V_esc, 2) - 2*r->G * initial_mass*(1/rtot - 1/params->separation_distance));

    if (big_frags == 1){  //assign radii, positions and velocities to second largest remnant
        struct reb_particle Slr1 = {0};
        Slr1.m = params->Mslr;
        Slr1.x = com.x - params->separation_distance*cos(theta_inc*(0))*unit_vix - params->separation_distance*sin(theta_inc*(0))*ox;
        Slr1.y = com.y - params->separation_distance*cos(theta_inc*(0))*unit_viy - params->separation_distance*sin(theta_inc*(0))*oy;
        Slr1.z = com.z - params->separation_distance*cos(theta_inc*(0))*unit_viz - params->separation_distance*sin(theta_inc*(0))*oz;

        Slr1.vx = com.vx - fragment_velocity*cos(theta_inc*(0))*unit_vix - fragment_velocity*sin(theta_inc*(0))*ox;
        Slr1.vy = com.vy - fragment_velocity*cos(theta_inc*(0))*unit_viy - fragment_velocity*sin(theta_inc*(0))*oy;
        Slr1.vz = com.vz - fragment_velocity*cos(theta_inc*(0))*unit_viz - fragment_velocity*sin(theta_inc*(0))*oz;
        Slr1.r = pow((3*Slr1.m)/(4*M_PI*params->rho), 1./3.);
        sprintf(hash,"FRAG%d", tot_no_frags+1);
        Slr1.hash = reb_hash(hash);
        Slr1.lastcollision = r->t;
        
        mxsum[0] = mxsum[0] + Slr1.m*Slr1.x;
        mxsum[1] = mxsum[1] + Slr1.m*Slr1.y;    
        mxsum[2] = mxsum[2] + Slr1.m*Slr1.z;

        mvsum[0] = mvsum[0] + Slr1.m*Slr1.vx;
        mvsum[1] = mvsum[1] + Slr1.m*Slr1.vy;   
        mvsum[2] = mvsum[2] + Slr1.m*Slr1.vz;

        reb_add(r, Slr1);
        Slr1.lastcollision = r->t;
    }

    int new_beginning_frag_index = tot_no_frags+big_frags+1;
    for (int i=(new_beginning_frag_index); i<(new_beginning_frag_index+no_frags); i++){          //add fragments
        struct reb_particle fragment = {0};
        if (frag_mass >= min_frag_mass){
        fragment.m = frag_mass;                  

        fragment.x = com.x - params->separation_distance*cos(theta_inc*(i))*unit_vix - params->separation_distance*sin(theta_inc*(i))*ox;
        fragment.y = com.y - params->separation_distance*cos(theta_inc*(i))*unit_viy - params->separation_distance*sin(theta_inc*(i))*oy;
        fragment.z = com.z - params->separation_distance*cos(theta_inc*(i))*unit_viz - params->separation_distance*sin(theta_inc*(i))*oz;

        fragment.vx = com.vx - fragment_velocity*cos(theta_inc*(i))*unit_vix - fragment_velocity*sin(theta_inc*(i))*ox;
        fragment.vy = com.vy - fragment_velocity*cos(theta_inc*(i))*unit_viy - fragment_velocity*sin(theta_inc*(i))*oy;
        fragment.vz = com.vz - fragment_velocity*cos(theta_inc*(i))*unit_viz - fragment_velocity*sin(theta_inc*(i))*oz;
        fragment.r = pow((3*fragment.m)/(4*M_PI*params->rho), 1./3.);
        sprintf(hash, "FRAG%d", i);
        

        fragment.hash = reb_hash(hash);
        fragment.lastcollision = r->t;
        mxsum[0] = mxsum[0] + fragment.m*fragment.x;
        mxsum[1] = mxsum[1] + fragment.m*fragment.y;    
        mxsum[2] = mxsum[2] + fragment.m*fragment.z;

        mvsum[0] = mvsum[0] + fragment.m*fragment.x;
        mvsum[1] = mvsum[1] + fragment.m*fragment.y;    
        mvsum[2] = mvsum[2] + fragment.m*fragment.z;

        reb_add(r, fragment); 
        fragment.lastcollision = r->t;
                                    }
                                }
    tot_no_frags += big_frags+no_frags;
    //Ensure momentum are conserved

    
    double xoff[3] = {com.x - mxsum[0]/initial_mass, com.y - mxsum[1]/initial_mass, com.z - mxsum[2]/initial_mass};
    double voff[3] = {com.vx - mvsum[0]/initial_mass, com.vy - mvsum[1]/initial_mass, com.vz - mvsum[2]/initial_mass};


    target -> x +=  xoff[0]; 
    target -> y += xoff[1]; 
    target -> z += xoff[2]; 
    target -> vx += voff[0]; 
    target -> vy += voff[1]; 
    target -> vz += voff[2]; 
    target->r = pow((3*target->m)/(4*M_PI*params->rho), 1./3.);

    for (int i=(tot_no_frags-new_bodies)+1; i<(tot_no_frags+1); i++){ 
        char frag[10];
        sprintf(frag, "FRAG%d", i);
        reb_get_particle_by_hash(r, reb_hash(frag))->x += xoff[0];
        reb_get_particle_by_hash(r, reb_hash(frag))->y += xoff[1];
        reb_get_particle_by_hash(r, reb_hash(frag))->z += xoff[2];

        reb_get_particle_by_hash(r, reb_hash(frag))->vx += voff[0];
        reb_get_particle_by_hash(r, reb_hash(frag))->vy += voff[1];
        reb_get_particle_by_hash(r, reb_hash(frag))->vz += voff[2];

    }

    return;
}


void elastic_bounce(struct reb_simulation* const r, struct reb_collision c, struct collision_params *params){
    struct reb_particle* particles = r->particles;

    int i = c.p1;   // i is target particle, j is projectile
    int j = c.p2; 

    if (particles[i].m<particles[j].m){        //particle with greater mass is the target; i = target, j = projectile
        i = c.p2;
        j = c.p1;
                }
    struct reb_particle com = reb_get_com_of_pair(particles[i], particles[j]);
    double msum, r_diff,p,vr,vtheta, vphi,temp, vrelx,vrely, vrelz;

    msum = particles[i].m + particles[j].m;
    //convert into spherical coordinates
    r_diff = sqrt(pow(params->xrel, 2));
    p = sqrt(pow(params->dx,2) + pow(params->dy,2));
    vr = get_dot(params->dx, params->dy, params->dz, params->Vix, params->Viy, params->Viz)/r_diff;

    vtheta = p*params->dz/r_diff - params->dz*(params->dx*params->Vix + params->dy*params->Viy)/(r_diff*p);
    vphi = (params->dx*params->Viy - params->dy*params->Vix)/p;

    //convert back into cartesian coordinates using -vr
    temp = p*-vr - params->dz*vtheta; 

    vrelx = params->dx*temp/(r_diff*p) - params->dy*vphi/p;
    vrely = params->dy*temp/(r_diff*p) + params->dx*vphi/p;
    vrelz = params->dz*vr + p*vtheta/r_diff;

    particles[i].x = com.x + particles[j].m*params->dx/msum;
    particles[i].y = com.y + particles[j].m*params->dy/msum;
    particles[i].z = com.z + particles[j].m*params->dz/msum;

    particles[j].x = com.x - particles[i].m*params->dx/msum;
    particles[j].y = com.y - particles[i].m*params->dy/msum;
    particles[j].z = com.z - particles[i].m*params->dz/msum;

    particles[i].vx = com.vx + particles[j].m*vrelx/msum;
    particles[i].vy = com.vy + particles[j].m*vrely/msum;
    particles[i].vz = com.vz + particles[j].m*vrelz/msum;

    particles[j].vx = com.vx - particles[i].m*vrelx/msum;
    particles[j].vy = com.vy - particles[i].m*vrely/msum;
    particles[j].vz = com.vz - particles[i].m*vrelz/msum;

    particles[i].x = com.x + particles[j].m*params->dx/msum + particles[i].vx*r->dt;
    particles[i].y = com.y + particles[j].m*params->dy/msum + particles[i].vy*r->dt;
    particles[i].z = com.z + particles[j].m*params->dz/msum + particles[i].vz*r->dt;

    particles[j].x = com.x - particles[i].m*params->dx/msum - particles[j].vx*r->dt;
    particles[j].y = com.y - particles[i].m*params->dy/msum - particles[j].vy*r->dt;
    particles[j].z = com.z - particles[i].m*params->dz/msum - particles[j].vz*r->dt;

    return;
}



int hit_and_run(struct reb_simulation* const r, struct reb_collision c, struct collision_params *params){  //also includes partial accretion.  Mlr = M_target.  Projectile is erroded.
        struct reb_particle* target = &(r->particles[params->target]);
        struct reb_particle* projectile = &(r->particles[params->projectile]);


        int swap = 2;
        int i = c.p1;
        int j = c.p2;   //make sure projectile is the particle being removed
        struct reb_particle* pi = &(r->particles[i]);
        struct reb_particle* pj = &(r->particles[j]);
        if (pi->m < pj->m){
            swap = 1;
                }

        double phi = 2*acos((params->l-projectile->r)/projectile->r);
        double A_interact = pow(projectile->r, 2)*((M_PI-(phi-sin(phi))/2.));  //Leinhardt Eq. 46;
        double L_interact = 2.*pow(pow(target->r,2)-(pow(target->r-params->l/2.,2)), .5);   //Leinhardt Eq. 47
        double beta = (A_interact*L_interact)/target->m;  //Chambers Eq. 11
        double Rc1 = pow(3./(4.*M_PI*params->rho1)*(beta*target->m + projectile->m), 1./3.); 
        double Q0 = .8*params->cstar*M_PI*params->rho1*r->G*pow(Rc1, 2);
        double gamma = (beta*target->m)/projectile->m;
        double Q_star = (pow(1+gamma, 2)/4*gamma)* Q0;

        double mu = (beta*target->m*projectile->m)/(beta*target->m+projectile->m);  //Chambers Eq. 13
        double Q = .5*(mu*pow(params->Vi,2))/(beta*target->m+projectile->m); //Chambers Eq. 12


        double c1 = 2.43;
        double c2 = -0.0408;
        double c3 = 1.86;
        double c4 = 1.08;

        double big_gamma = pow((1-gamma)/(1+gamma), 2);
        double theta = acos(params->vrel/params->xrel);
        double v_crit = params->V_esc*(c1*big_gamma*(pow(1-sin(theta), 5./2.)) + c2*big_gamma + c3*pow(1-sin(theta), 5./2.) + c4);
        double targ_m = target->m;
        double imp_m = projectile->m;


        if (params->Vi <= v_crit || targ_m + imp_m - params->Mlr < min_frag_mass){             //if impact velocity is low, the hit-and-run results in a merger.
            reb_collision_resolve_merge(r,c);
            return swap;
        }


        if (params->Mlr < targ_m){  //erosion of target with fragments            
            double new_Mlr = MAX(params->Mlr, min_frag_mass);
            params->Mlr = new_Mlr;
            params->Mslr = 0;
            add_fragments(r,c,params);
            return swap;
        }

        double Mlr_dag = (beta*target->m + projectile->m)/10 * pow(Q/(1.8*Q_star), -1.5);

        if (Q < 1.8*Q_star){
            Mlr_dag = (beta*targ_m + imp_m)*(1 - Q/ (2*Q_star));
        }

        if (Mlr_dag < min_frag_mass || (imp_m - Mlr_dag) < min_frag_mass){ 
            elastic_bounce(r,c, params);
            swap = 0;
            return swap;} //projectile is not erroded, both particles remain in sim

        else{
        params->Mslr = Mlr_dag;
        add_fragments(r,c,params);
        }

        return swap;
    }

void init_collision_params(struct collision_params* params){
    params->target=0;
    params->projectile=0;
    params->dx=0;
    params->dy=0;
    params->dz=0;
    params->b=0;
    params->Vix=0;
    params->Viy=0;
    params->Viz=0;
    params->Vi=0;
    params->l=0;
    params->rho1=0;
    params->rho=0;
    params->cstar=0;
    params->mubar=0;
    params->mu=0;
    params->QR=0;
    params->QpRD=0;
    params->V_esc=0;
    params->separation_distance=0;
    params->Mlr=0;
    params->Mslr=0;
    params->Q=0;
    params->Mlr_dag=0;
    params->Q_star=0;
    params->vrel=0;
    params->xrel=0;
}

struct collision_params* create_collision_params(){
    struct collision_params* params = calloc(1,sizeof(struct collision_params));
    init_collision_params(params);
    return params;
}

int get_collision_type(struct reb_simulation* const r, struct reb_collision c){
    struct reb_particle* particles = r->particles;
    struct collision_params* params = create_collision_params();
    const double G = r->G;
    int swap = 2;
    int i = c.p1;   // i is target particle, j is projectile
    int j = c.p2; 

    if (r->particles[i].m < r->particles[j].m){        //unless swap is redfined as 0, projectile is going to be removed.
        swap =1;
        i = c.p2;
        j = c.p1;
                }


    if (r->particles[i].lastcollision==r->t || r->particles[j].lastcollision==r->t) return 0;


    double imp_r = particles[j].r;
    double targ_r = particles[i].r;
    double R_tot = imp_r + targ_r;
    double imp_m = particles[j].m;
    double targ_m = particles[i].m;
    double M_tot = imp_m + targ_m;

    double Mlr,dx,dy,dz,Vix,Viy,Viz;
    double x2rel, xrel, v2rel, v2imp, Vi;
    double hx,hy,hz,h2,b;
    make_vector(particles[i].x, particles[i].y, particles[i].z, particles[j].x, particles[j].y, particles[j].z, &dx,&dy,&dz);  //find relative coordinates dx, dy,dz
    x2rel = get_dot(dx,dy,dz,dx,dy,dz);  
    make_vector(particles[i].vx, particles[i].vy, particles[i].vz, particles[j].vx, particles[j].vy, particles[j].vz, &Vix,&Viy,&Viz);  //find relative velocity
    v2rel = get_dot(Vix,Viy,Viz,Vix,Viy,Viz);

    xrel = sqrt(x2rel);  //distance between the centers of the projectile and target


    hx = (dy*Viz - dz*Viy);                     //angular momentum vector xrel X Vrel
    hy = (dz*Vix - dx*Viz);
    hz = (dx*Viy - dy*Vix);

    h2 = get_dot(hx,hy,hz,hx,hy,hz);

    v2imp = v2rel + 2*G*M_tot*(1./R_tot - 1./xrel); //impact velocity with gravitational focusing at time of detected collision

    if (1./R_tot - 1./xrel < 0){v2imp = v2rel;}  //if collision is detected after physical contact
    
    Vi = sqrt(v2imp);  //magnitude of impact velocity vector

    b = sqrt(h2/v2imp);  //impact parameter
    //the impact paramater is the distance between the two centers at the time of the collision, projected perpendicular to the impact velocity vector. (Wallace et al 2018)
    
    //Stewart & Leinhardt 2012 parameters
    double l = R_tot-b;  //Leinhardt Eq. 7, the projected length of the projectile overlapping the target
    l = MIN(l, 2*imp_r);
    double alpha = (pow(l,2)*(3*imp_r-l))/(4*pow(imp_r, 3)); //Leinhardt Eq. 11, interacting mass fraction
    alpha = MIN(1, alpha);
    double Q = .5*v2imp*targ_m*imp_m/pow(M_tot,2);  //specific energy per unit mass

    double v2esc = 2*G*M_tot/R_tot; //mutal escape velocity as defined in Wallace et al 2018 and Chambers 2013
    double V_esc = pow(v2esc, .5); //mutal escape velocity as defined in Wallace et al 2018 and Chambers 2013
    double mu = (targ_m*imp_m)/M_tot;  //Chambers Eq. 2, reduced mass
    double alphamu = (alpha*targ_m*imp_m)/(alpha*imp_m + targ_m);  //Leinhardt Eq. 12, reduced interacting mass for fraction alpha.
    double gamma = imp_m/targ_m;  //Chambers Eq. 6
    if (gamma > 1){gamma = targ_m/imp_m;}

    const double cstar = 1.8;      //user defined variable, default taken from paper
    const double mubar = 0.36;    // user defined variable, default taken from paper
    double rho1 = 1000;         //constant density in kg/cm^3, DO NOT CHANGE
    if (r->G == 1){rho1 = 1.684e12;}
    


    double Rc1 = pow((M_tot*3)/(4*M_PI*rho1), 1./3.);  //Chambers Eq. 4, combined radius of target and projectile with constant density
    double Q0 = .8*cstar*M_PI*rho1*G*pow(Rc1,2);  //Chambers Eq. 3, critical value of impact energy for head-on collisions
    double Q_star = pow(mu/alphamu, 1.5)*(pow(1+gamma, 2)/ 4*gamma)*Q0;  //Chambers Eq. 5, critical value for oblique or different mass collisons.  
    if (alpha == 0.0){Q_star = 6364136223846793005.0;}
    if (b == 0 && imp_m == targ_m){
        Q_star = Q0;
    }
    
    double qratio = Q/Q_star;
    if (qratio < 1.8){
        Mlr = M_tot*(1.0-.5*qratio);
    }
    else{
        Mlr = .1*M_tot*pow(qratio/1.8, -1.5);  //Chambers Eq.8
    }


    double separation_distance = 4 * (targ_r + imp_r);  //seperation distance of fragments.  Should be userdefined but this is what chambers uses
///POPULATE STRUCT OBJECTS
    params->target = i;
    params->projectile =j;
    params->dx = dx;
    params->dy = dy;
    params->dz = dz;
    params->b = b;
    params->Vix = Vix;
    params->Viy = Viy;
    params->Viz = Viy;
    params->Vi = Vi;
    params->l = l;
    params->rho1 = rho1;
    params->cstar = cstar;
    params->mubar = mubar;
    params->mu = mu;
    params->Q = Q;
    params->rho = 5.05e6;  //3 g/cm^3 in solar mass/AU^3
    params->separation_distance = separation_distance;
    params->V_esc = V_esc;
    params->vrel = sqrt(v2rel);
    params->Mslr = 0;
    params->xrel = xrel;
    params->Mlr = Mlr;


    if (Vi <= V_esc){
        reb_collision_resolve_merge(r,c);
        return swap;
    }
    else{
        if (b < targ_r){   //non-grazing regime  
            if (M_tot - params->Mlr < min_frag_mass){
                reb_collision_resolve_merge(r,c);
                return swap;
            }
            else {
            params->Mlr = MAX(Mlr, min_frag_mass);
            add_fragments(r,c,params);
            }
        }
        if (b >= targ_r){  //Chambers Eq. 9, hit and run criteria, grazing regime
            swap = hit_and_run(r,c,params);
        }
        }
    return swap;    
}



int reb_collision_resolve_fragment(struct reb_simulation* const r, struct reb_collision c){
    if (r->particles[c.p1].lastcollision==r->t || r->particles[c.p2].lastcollision==r->t) return 0;
    int i = c.p1;
    int j = c.p2; 

    if (i<j) return 0;      //only return one collision callback

    struct reb_particle* pj = &(r->particles[j]);
    struct reb_particle* pi = &(r->particles[i]);

    if (pj->m >= mass_limit){
        reb_collision_resolve_merge(r,c);
        return 1;
    }

    if (pi->m >= mass_limit){
        reb_collision_resolve_merge(r,c);
        return 2;
    }  
}
