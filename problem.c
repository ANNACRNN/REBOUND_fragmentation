#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "/home/acc298/rebound/src/rebound.h"



double min_frag_mass = 1.4e-8;
int tot_no_frags = 0;  //if restarting a simulation this needs to be changed to the last number of frags in the simulation, otherwise new fragments added will rewrite exisiting frags


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
    double cstar;
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
    int collision_type;
    int no_frags;
}; 


void make_vector(double x1, double y1, double z1, double x2, double y2, double z2, double *x, double*y, double*z){   //Galilean transform
    *x = x1-x2;
    *y = y1-y2;
    *z = z1-z2;
}

double get_dot(double x1, double y1, double z1, double x2, double y2, double z2){ 
    return (x1*x2)+(y1*y2)+(z1*z2);
        }  //return dot product of two vectors

double get_mag(double x, double y, double z){
    return sqrt(pow(x,2)+pow(y,2)+pow(z,2));
        }   //return magnitude of vector
double get_radii(double m, double rho){
    return pow((3*m)/(4*M_PI*rho),1./3.);
}

void add_fragments(struct reb_simulation* const r, struct reb_collision c, struct collision_params *params){
    struct reb_particle* target = &(r->particles[params->target]);
    struct reb_particle* projectile = &(r->particles[params->projectile]);
    struct reb_particle com = reb_particle_com_of_pair(*target, *projectile);
    double initial_mass = target -> m + projectile -> m;
    double remaining_mass = initial_mass - params->Mlr;
    double rho = target->m/(4./3*M_PI*pow(target ->r, 3));
    double rtot = target -> r + projectile -> r;

    int big_frags = 0;
    if (params->Mslr > 0){
        remaining_mass = remaining_mass -  params->Mslr;
        big_frags = 1;
    }

    int no_frags = remaining_mass/min_frag_mass;  //fragments are broken up into equal sizes
    double frag_mass = remaining_mass/no_frags;


    int new_bodies = no_frags + big_frags;
    params->no_frags = new_bodies;
    

    char hash[10];
    double mxsum[3] = {0,0,0};
    double mvsum[3] = {0,0,0};
    //target gets mass of Mlr and is assigned COM position and velocity;
    target -> last_collision = r->t;
    target -> m = params->Mlr;
    target -> r = get_radii(params->Mlr, rho);
    target->x = com.x;
    target->y = com.y;
    target->z = com.z;

    target->vx = com.vx;
    target->vy = com.vy;
    target->vz = com.vz;
    
    if (no_frags == 1 && params->Mlr <= frag_mass){
        target->m = frag_mass;
        frag_mass = params->Mlr;
    }

    mxsum[0] = mxsum[0] + target->m*target->x;
    mxsum[1] = mxsum[1] + target->m*target->y;  
    mxsum[2] = mxsum[2] + target->m*target->z;

    mvsum[0] = mvsum[0] + target->m*target->vx;
    mvsum[1] = mvsum[1] + target->m*target->vy; 
    mvsum[2] = mvsum[2] + target->m*target->vz;         
    
    double theta_inc = (2.*M_PI)/new_bodies;
    

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

    double fragment_velocity =sqrt(1.1*pow(params->V_esc,2) - 2*r->G*initial_mass*(1./rtot - 1./params->separation_distance));

    if (big_frags == 1){  //assign radii, positions and velocities to second largest remnant, theta=0
        struct reb_particle Slr1 = {0};
        Slr1.m = params->Mslr;
        Slr1.x = com.x + params->separation_distance*unit_vix;
        Slr1.y = com.y + params->separation_distance*unit_viy;
        Slr1.z = com.z + params->separation_distance*unit_viz;

        Slr1.vx = com.vx + fragment_velocity*unit_vix;
        Slr1.vy = com.vy + fragment_velocity*unit_viy;
        Slr1.vz = com.vz + fragment_velocity*unit_viz;

        Slr1.r = get_radii(Slr1.m, rho);
        sprintf(hash,"FRAG%d", tot_no_frags+1);
        Slr1.hash = reb_hash(hash);
        printf("%s hash, mass:      %u %e\n", hash, Slr1.hash, Slr1.m);
        mxsum[0] += Slr1.m*Slr1.x;
        mxsum[1] += Slr1.m*Slr1.y;    
        mxsum[2] += Slr1.m*Slr1.z;

        mvsum[0] += Slr1.m*Slr1.vx;
        mvsum[1] += Slr1.m*Slr1.vy;   
        mvsum[2] += Slr1.m*Slr1.vz;
        Slr1.last_collision = r->t;
        reb_simulation_add(r, Slr1);
    }



    int new_beginning_frag_index = tot_no_frags+big_frags+1;
    for (int i=(new_beginning_frag_index); i<(new_beginning_frag_index+no_frags); i++){          //add fragments
        struct reb_particle fragment = {0};
        int j = i - new_beginning_frag_index+1;
        fragment.m = frag_mass;                  
        fragment.x = com.x + params->separation_distance*(cos(theta_inc*j)*unit_vix + sin(theta_inc*j)*ox);
        fragment.y = com.y + params->separation_distance*(cos(theta_inc*j)*unit_viy + sin(theta_inc*j)*oy);
        fragment.z = com.z + params->separation_distance*(cos(theta_inc*j)*unit_viz + sin(theta_inc*j)*oz);
        fragment.vx = com.vx + fragment_velocity*(cos(theta_inc*j)*unit_vix + sin(theta_inc*j)*ox);
        fragment.vy = com.vy + fragment_velocity*(cos(theta_inc*j)*unit_viy + sin(theta_inc*j)*oy);
        fragment.vz = com.vz + fragment_velocity*(cos(theta_inc*j)*unit_viz + sin(theta_inc*j)*oz);

        fragment.r = get_radii(frag_mass, rho);
        fragment.last_collision = r->t;
        sprintf(hash, "FRAG%d", i);
        fragment.hash = reb_hash(hash);
        printf("%s hash, mass:      %u %e\n", hash, fragment.hash, fragment.m);
        mxsum[0] +=fragment.m*fragment.x;
        mxsum[1] += fragment.m*fragment.y;    
        mxsum[2] += fragment.m*fragment.z;

        mvsum[0] += fragment.m*fragment.vx;
        mvsum[1] += fragment.m*fragment.vy;    
        mvsum[2] += fragment.m*fragment.vz;

        reb_simulation_add(r, fragment); 
                                }
    tot_no_frags += big_frags+no_frags;


    //Ensure momentum is conserved


    
    double xoff[3] = {com.x - mxsum[0]/initial_mass, com.y - mxsum[1]/initial_mass, com.z - mxsum[2]/initial_mass};
    double voff[3] = {com.vx - mvsum[0]/initial_mass, com.vy - mvsum[1]/initial_mass, com.vz - mvsum[2]/initial_mass};


    target -> x +=  xoff[0]*target->m/initial_mass; 
    target -> y += xoff[1]*target->m/initial_mass; 
    target -> z += xoff[2]*target->m/initial_mass; 
    target -> vx += voff[0]*target->m/initial_mass; 
    target -> vy += voff[1]*target->m/initial_mass; 
    target -> vz += voff[2]*target->m/initial_mass; 

    for (int i=(tot_no_frags-new_bodies)+1; i<(tot_no_frags+1); i++){ 
        char frag[10];
        sprintf(frag, "FRAG%d", i);
        double mass_fraction = reb_simulation_particle_by_hash(r, reb_hash(frag))->m/initial_mass;
        reb_simulation_particle_by_hash(r, reb_hash(frag))->x += xoff[0]*mass_fraction;
        reb_simulation_particle_by_hash(r, reb_hash(frag))->y += xoff[1]*mass_fraction;
        reb_simulation_particle_by_hash(r, reb_hash(frag))->z += xoff[2]*mass_fraction;

        reb_simulation_particle_by_hash(r, reb_hash(frag))->vx += voff[0]*mass_fraction;
        reb_simulation_particle_by_hash(r, reb_hash(frag))->vy += voff[1]*mass_fraction;
        reb_simulation_particle_by_hash(r, reb_hash(frag))->vz += voff[2]*mass_fraction;
    }

    return;
}


void merge(struct reb_simulation* const r, struct reb_collision c, struct collision_params *params){
    struct reb_particle* pi = &(r->particles[params->target]);
    struct reb_particle* pj = &(r->particles[params->projectile]);

    double invmass = 1.0/(pi->m + pj->m);
    double targ_rho = pi->m/(4./3*M_PI*pow(pi->r,3));  //new body recieves density of the target
    // Merge by conserving mass, volume and momentum
    pi->vx = (pi->vx*pi->m + pj->vx*pj->m)*invmass;
    pi->vy = (pi->vy*pi->m + pj->vy*pj->m)*invmass;
    pi->vz = (pi->vz*pi->m + pj->vz*pj->m)*invmass;
    pi->x  = (pi->x*pi->m + pj->x*pj->m)*invmass;
    pi->y  = (pi->y*pi->m + pj->y*pj->m)*invmass;
    pi->z  = (pi->z*pi->m + pj->z*pj->m)*invmass;
    pi->m  = pi->m + pj->m;
    pi->r  = pow((3*pi->m)/(4*M_PI*targ_rho),1./3.);
    pi->last_collision = r->t;


    return; // 
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

        double targ_m = target->m;
        double imp_m = projectile->m;
        double zeta = pow((targ_m - imp_m)/(targ_m + imp_m),2);
        double fac = pow(1-params->b/(target->r + projectile->r),2.5);
        double v_crit = params->V_esc*(c1*zeta*fac + c2*zeta +c3*fac + c4);

        if (params->Vi <= v_crit){             //if impact velocity is low, the hit-and-run results in a merger.
            printf("GRAZE AND MERGE\n");  
            params->collision_type = 1;          
            merge(r,c,params);
            return swap;
        }

        else{ //vi>v_crit
            if (params->Mlr<targ_m){ //Target is being eroded, projectile should also fragment
                if (targ_m+imp_m <= 2*min_frag_mass){ //not enough mass to produce new fragments
                    printf("ELASTIC BOUNCE\n");
                    params->collision_type=0;
                    reb_collision_resolve_hardsphere(r,c);
                    swap = 0;
                                                       }
                else{
                    params->Mlr = MAX(params->Mlr, min_frag_mass);
                    printf("GRAZING PARTIAL EROSION\n");
                    params->collision_type = 3;
                    add_fragments(r,c,params);
                    }
                                    }
            else{ //Mlr > Mt, either a hit and run or an elastic bounce
                double Mlr_dag = (beta*target->m + projectile->m)/10 * pow(Q/(1.8*Q_star), -1.5);
                if (Q < 1.8*Q_star){
                    Mlr_dag = (beta*targ_m + imp_m)*(1 - Q/ (2*Q_star));
                }

            double projectile_mass_accreted = params->Mlr - targ_m;
            double new_projectile_mass = projectile->m - projectile_mass_accreted;
            Mlr_dag = MAX(Mlr_dag, min_frag_mass);
            if (new_projectile_mass-Mlr_dag < min_frag_mass){
                    printf("ELASTIC BOUNCE\n");
                    params->collision_type=0;
                    reb_collision_resolve_hardsphere(r,c);
                    swap = 0;
                                                            }
                else{
                    params->Mslr = Mlr_dag;
                    printf("HIT AND RUN\n");
                    params->collision_type = 2;
                    add_fragments(r,c,params);
                    }
                }
        return swap;
        }
    }

void print_collision_array(struct reb_simulation* const r, struct reb_collision c, struct collision_params *params){  
//0=elastic bounce, 1=merger, 2=partial accretion, 3=partial erosion, 4=supercat
    FILE* of = fopen("collision_report.txt","a+");
    fprintf(of, "%e\t", r->t);     
    fprintf(of, "%d\t", params->collision_type);                  
    fprintf(of, "%u\t", (r->particles[params->target].hash));  
    fprintf(of, "%e\t", (r->particles[params->target].m));
    fprintf(of, "%u\t", (r->particles[params->projectile].hash)); 
    for(int i=(r->N - params->no_frags);i<r->N;i++){        //assuming Fragments are added to end of particle array
        fprintf(of, "%u\t", (r->particles[i].hash));  
        fprintf(of, "%e\t", (r->particles[i].m)); 
    }
    fprintf(of, "\n");   
    fclose(of);                        // close file

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
    params->cstar=0;
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
    params->collision_type=0;
    params->no_frags = 0;
}

struct collision_params* create_collision_params(){
    struct collision_params* params = calloc(1,sizeof(struct reb_simulation));
    init_collision_params(params);
    return params;
}


int reb_collision_resolve_fragment(struct reb_simulation* const r, struct reb_collision c){
    if (r->particles[c.p1].last_collision==r->t || r->particles[c.p2].last_collision==r->t) return 0;
    int i = c.p1;
    int j = c.p2; 
    if (i<j) return 0;      //only return one collision callback

    int swap = 2;
    if (r->particles[i].m < r->particles[j].m){        //unless swap is redfined as 0, projectile is going to be removed.
        swap =1;
        i = c.p2;
        j = c.p1;
                }

    struct reb_particle* particles = r->particles;
    struct collision_params* params = create_collision_params();

    double imp_r = particles[j].r;
    double targ_r = particles[i].r;
    double R_tot = imp_r + targ_r;

    double imp_m = particles[j].m;
    double targ_m = particles[i].m;

    printf("TIME OF COLLISION: %e\n", r->t);
    printf("Target hash, mass = %u %e\n", particles[i].hash, targ_m);
    printf("Projectile hash, mass = %u %e\n", particles[j].hash, imp_m);

    double M_tot = imp_m + targ_m;
    double G = r->G;
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
    b = sqrt(h2/v2imp);  //impact parameter, b=R_tot*sin(theta)
    if (b != b){
        printf("NAN b \n");
        exit(0);}
    //Stewart & Leinhardt 2012 parameters
    double mu = (targ_m*imp_m)/M_tot;  //Chambers Eq. 2, reduced mass
    double l = R_tot-b;  //Leinhardt Eq. 7, the projected length of the projectile overlapping the target
    l = MIN(l, 2*imp_r);
    double alpha = (pow(l,2)*(3*imp_r-l))/(4*pow(imp_r, 3)); //Leinhardt Eq. 11, interacting mass fraction
    alpha = MIN(1., alpha);
    double Q = .5*v2imp*targ_m*imp_m/pow(M_tot,2);  //specific energy per unit mass
    double V_esc = pow(2.*G*M_tot/R_tot, .5); //mutal escape velocity as defined in Wallace et al 2018 and Chambers 2013
    double alphamu = (alpha*targ_m*imp_m)/(alpha*imp_m + targ_m);  //Leinhardt Eq. 12, reduced interacting mass for fraction alpha.
    double gamma = imp_m/targ_m;  //Chambers Eq. 6

    const double cstar = 1.8;      //may be a user defined variable, default taken from paper

    double rho1;         //constant density

    if (G==6.674e-8){rho1 =1;} //CGS
    if (G==6.674e-11){rho1 =1000;} //SI
    if (G==39.476926421373 || G==1){rho1 = 1.684e6;}  //Msun/AU^3
    double Rc1 = pow((M_tot*3)/(4.*M_PI*rho1), 1./3.);  //Chambers Eq. 4, combined radius of target and projectile with constant density
    double Q0 = .8*cstar*M_PI*rho1*G*pow(Rc1,2);  //Chambers Eq. 3, critical value of impact energy for head-on collisions
    double Q_star = pow(mu/alphamu, 1.5)*(pow(1+gamma, 2)/ (4*gamma))*Q0;  //Chambers Eq. 5, critical value for oblique or different mass collisons.  
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

    double separation_distance = 4 * R_tot;  //seperation distance of fragments.  Should be userdefined but this is what chambers uses
///POPULATE STRUCT OBJECTS
    params->target = i;
    params->projectile =j;
    params->dx = dx;
    params->dy = dy;
    params->dz = dz;
    params->b = b;
    params->Vix = Vix;
    params->Viy = Viy;
    params->Viz = Viz;
    params->Vi = Vi;
    params->l = l;
    params->rho1 = rho1;
    params->cstar = cstar;
    params->mu = mu;
    params->Q = Q;
    params->separation_distance = separation_distance;
    params->V_esc = V_esc;
    params->vrel = sqrt(v2rel);
    params->Mslr = 0;
    params->xrel = xrel;
    params->Mlr = Mlr;


    printf("Mp/Mt:    %0.4f\n", imp_m/targ_m);
    printf("Mlr/Mt:    %0.4f\n", Mlr/targ_m);
    printf("Mlr/Mtot:    %0.4f\n", Mlr/M_tot);
    printf("b/Rtarg:     %0.4f\n", b/targ_r);
    printf("Vimp/Vesc:     %0.4f\n",  Vi/V_esc);
    printf("Q/ Qstar:     %0.4f\n", Q/Q_star);
    printf("COLLISION TYPE: ");

    if (Vi <= V_esc){
        params->collision_type = 1;
        printf("SIMPLY MERGED\n");
        merge(r,c, params);
                    }
    else{  //Vi > V_esc
        if (b<targ_r){ //non-grazing regime
            if (M_tot - params->Mlr < min_frag_mass){
                params->collision_type = 1;
                printf("EFFECTIVELY MERGED\n");
                merge(r,c,params);
                print_collision_array(r,c,params);  
                                                     }
            else{ // M_tot - params->Mlr >= min_frag_mass; fragments will be produced unless it is a graze and merge or elastic bounce 
                if (params->Mlr < targ_m){
                    if (params->Mlr <= 0.1*targ_m){
                        printf("SUPER-CATASTROPHIC\n");
                        params->collision_type = 4;
                        params->Mlr = MAX(Mlr, min_frag_mass);
                        add_fragments(r,c,params);
                                                   }
                    else{
                        printf("PARTIAL EROSION\n");
                        params->collision_type = 3;
                        params->Mlr = MAX(Mlr, min_frag_mass);
                        add_fragments(r,c,params);
                        }
                                          }
                                    
                else{  //(params->Mlr >= targ_m)
                            printf("PARTIAL ACCRETION\n");
                            params->collision_type = 2;
                            add_fragments(r,c,params);
                    }
                }
                     }
        else{ // b > b_crit, grazing regime
            swap = hit_and_run(r,c,params); //swap gets redefined here as it may be set to 0 in the case of a bounce
            }
        }
                                

print_collision_array(r,c,params);  
return swap;
}


void heartbeat(struct reb_simulation* r);

double masses[288]={
2.791e-07, 2.791e-07, 2.791e-07, 2.791e-07, 2.791e-07, 2.791e-07, 2.791e-07, 2.791e-07, 2.791e-07, 2.791e-07, 2.791e-07, 2.791e-07, 2.791e-07, 2.791e-07, 2.791e-07, 2.791e-07, 2.791e-07, 2.791e-07, 2.791e-07, 2.791e-07, 2.791e-07, 2.791e-07, 2.791e-07, 2.791e-07, 2.791e-07, 2.791e-07, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08, 2.791e-08};
double semi_major[288]={
0.458968,0.665599,0.756319,0.844885,0.938353,1.036723,1.139997,1.248172,1.36125,1.479231,1.602114,1.729899,1.862587,2.000178,2.1425242,2.2898926,2.4421612,2.59933,2.761399,2.9283682,3.1002376,3.2770072,3.458677,3.645247,3.8367172,4.0330876,0.356999,0.433951,0.468586,0.495262,0.517556,0.536972,0.554315,0.570075,0.584581,0.598059,0.610679,0.622568,0.633825,0.64453,0.654746,0.664528,0.673919,0.682956,0.691672,0.700096,0.708406,0.716764,0.725172,0.733629,0.742135,0.75069,0.759293,0.767946,0.776648,0.785399,0.794199,0.803048,0.811946,0.820893,0.829889,0.838934,0.848028,0.857171,0.866363,0.875605,0.884895,0.894234,0.903622,0.91306,0.922546,0.932081,0.941665,0.951299,0.960981,0.970713,0.980493,0.990323,1.000201,1.010129,1.020105,1.030131,1.040205,1.050329,1.060501,1.070723,1.080994,1.091314,1.101682,1.1121,1.122567,1.133083,1.143647,1.154261,1.164924,1.175636,1.186397,1.197207,1.208066,1.218974,1.229931,1.240937,1.251992,1.263096,1.274249,1.285452,1.296703,1.308003,1.319352,1.330751,1.342198,1.353694,1.365239,1.376834,1.388477,1.40017,1.411911,1.423701,1.435541,1.447429,1.459367,1.471354,1.483389,1.495474,1.507607,1.51979,1.532022,1.544302,1.556632,1.569011,1.581439,1.593916,1.606441,1.619016,1.63164,1.644313,1.657035,1.669806,1.682626,1.695495,1.708413,1.72138,1.734396,1.747461,1.760575,1.773738,1.786951,1.800212,1.813522,1.826881,1.84029,1.853747,1.867253,1.880809,1.894413,1.908066,1.921769,1.93552,1.949321,1.96317,1.977069,1.991016,2.005013,2.019059,2.033153,2.047297,2.05317058,2.06741352,2.08170551,2.09604652,2.11043657,2.12487565,2.13936377,2.15390091,2.16848709,2.18312231,2.19780655,2.21253983,2.22732215,2.24215349,2.25703387,2.27196328,2.28694173,2.30196921,2.31704572,2.33217126,2.34734584,2.36256945,2.3778421,2.39316377,2.40853448,2.42395423,2.439423,2.45494081,2.47050766,2.48612353,2.50178844,2.51750238,2.53326536,2.54907737,2.56493841,2.58084848,2.59680759,2.61281573,2.6288729,2.64497911,2.66113435,2.67733862,2.69359193,2.70989427,2.72624564,2.74264604,2.75909548,2.77559395,2.79214146,2.808738,2.82538357,2.84207817,2.85882181,2.87561448,2.89245618,2.90934692,2.92628669,2.94327549,2.96031333,2.97740019,2.9945361,3.01172103,3.028955,3.046238,3.06357003,3.0809511,3.0983812,3.11586034,3.1333885,3.1509657,3.16859193,3.1862672,3.2039915,3.22176483,3.2395872,3.2574586,3.27537903,3.29334849,3.31136699,3.32943452,3.34755108,3.36571668,3.38393131,3.40219497,3.42050767,3.4388694,3.45728016,3.47573996,3.49424879,3.51280665,3.53141354,3.55006947,3.56877443,3.58752842,3.60633145,3.62518351,3.64408461,3.66303473,3.68203389,3.70108208,3.72017931,3.73932557,3.75852086,3.77776519,3.79705854,3.81640094,3.83579236,3.85523282,3.87472231,3.89426083,3.91384839,3.93348498,3.9531706,3.97290526,3.99268895,4.01252167,4.03240342,4.05233421,4.07231403,4.09234289};


int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();
    r->G = 39.476926421373;
    r->dt = 6./365.;
    r->integrator = REB_INTEGRATOR_MERCURIUS;
    r->heartbeat = heartbeat;
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_fragment;
    r->rand_seed = 1;

        //Assigning mass and number of planetary embryos and planetesimals
    double rho = 5.05e6; //3 g/cm^3


    struct reb_particle star = {0};
    star.m = 1.00;
    star.r = 0.1; 
    star.hash = 0; 
    reb_simulation_add(r, star);
    
    // Add planetary embryos
    
    for (int i=0; i<286; i++){
        double a = semi_major[i];      // semi major axis
        double m = masses[i];
        double inc = reb_random_uniform(r,0,0.0175);
        double ecc = reb_random_uniform(r,0,0.01);
        double omega = reb_random_uniform(r,0,2*M_PI);
        double Omega = reb_random_uniform(r,0,2*M_PI);
        double f = reb_random_uniform(r,0,2*M_PI);
        double hash = i + 1;
        //now build particle from orbit
        struct reb_particle emb = reb_particle_from_orbit(r->G, star, m, a, ecc, inc, Omega, omega, f);
        emb.r = get_radii(m, rho)*5;
        emb.hash = hash;
        
        reb_simulation_add(r, emb); 
    }
        double m,a,e,inc,Omega,omega,f; //Omega=longitude of ascending node, omega= argument of pericenter in RADIANS

    //Add Jupiter and Saturn
    struct reb_particle Jup = reb_particle_from_orbit(r->G, r->particles[0], m=9.543e-4, a=5.20349, e=0.048381, inc=0.365*(M_PI/180), Omega=0.0, omega=68.3155*(M_PI/180), f=227.0537*(M_PI/180));
    Jup.r = get_radii(Jup.m, rho); 
    Jup.hash = reb_hash("JUPITER");
    reb_simulation_add(r,Jup);

    struct reb_particle Sat = reb_particle_from_orbit(r->G, r->particles[0], m=0.0002857, a=9.54309, e=0.052519, inc=0.8892*(M_PI/180), Omega=M_PI, omega=324.5263*(M_PI/180),f=256.9188*(M_PI/180));
    Sat.r = get_radii(Sat.m, rho);
    Sat.hash = reb_hash("SATURN");
    reb_simulation_add(r,Sat);


    reb_simulation_move_to_com(r);  // This makes sure the planetary systems stays within the computational domain and doesn't drift.
    reb_simulation_save_to_file_interval(r,"simulationarchive.bin",1.e5);
    double run_time = 5.e6;
    reb_simulation_integrate(r, run_time);

}

void heartbeat(struct reb_simulation* r){
    if (reb_simulation_output_check(r, 1.e3)){  
        reb_simulation_output_timing(r, 0);
        printf("Walltime(s) = %f \n", r->walltime); 
    }
}
