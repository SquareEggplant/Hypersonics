//CYLINDRICAL WAVERIDER BASED OFF TOP SURFACE FUNCTION
//ALL ISSUES ARE FAULT OF ERIC JO
//Last updated 2025/10/16

#include <cmath>
#include <vector>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <sstream>

namespace fs = std::filesystem;

struct Vcone{
    double Vr;
    double Vtheta; 
};

struct vec{
    double x;
    double y;
    double z;
};

struct tri{
    int x;
    int y;
    int z;
    int comp;
};

struct limits{
    double x_edge;
    double length;
    double planar_area;
};

double TMS(double Vr, double Vtheta, double gamma, double Shock_Angle){
    //Euler Explicit Solution of the Taylor Maccoll Equation
    double gamma_factor = (gamma-1)/2;

    double numer = gamma_factor * (1 - (Vr*Vr) - (Vtheta*Vtheta)) * ((2*Vr) + (Vtheta/tan(Shock_Angle))) - (Vr*(Vtheta*Vtheta));
    double denom = (Vtheta*Vtheta) - gamma_factor * (1 - (Vr*Vr) - (Vtheta*Vtheta));

    double F2 = numer/denom;
    return F2;
}

void readinputs(std::istream& inp, std::string& line, std::string varname){
    const char delim = 58;
    while(line!=varname && !inp.eof()){
        std::getline(inp, line, delim);
    }
    std::getline(inp, line, delim);
}

Vcone ICs(double beta, double M1, double gamma){
    //Eric Jo's Initial Condition Function
    //Normal free-stream mach (M1) 
    double sinb = sin(beta);
    double Mn1 = M1*sinb;

    //Normal M2
    double Mn2 = sqrt( (1+0.5*(gamma-1)*Mn1*Mn1) / (gamma*Mn1*Mn1 - 0.5*(gamma-1)) );

    //Theta-Beta-Mach Wedge Estimation
    double tangent_theta = 2/tan(beta)*(((M1*M1)*(sinb*sinb) - 1)/((M1*M1)*(gamma+cos(2*beta)) + 2));
    double theta = atan(tangent_theta);

    //M2 
    double M2 = Mn2/sin(beta-theta);

    //Non-dimensional initial velocities behind the shock
    double V2 = 1/sqrt((2/(gamma-1))/M2/M2 + 1);
    double V_r0 = V2*cos(beta-theta);
    double V_theta0 = -V2*sin(beta-theta);

    Vcone VIC;
    VIC.Vr = V_r0;
    VIC.Vtheta = V_theta0;

    return VIC;
}

Vcone Interp(double theta, std::vector<Vcone> Vel, std::vector<double> thetavals){
    int ind;
    Vcone res;
    double d, dx, dyv, dyt;

    for (int i=0; i<thetavals.size();i++){
        if(thetavals[i]<theta){
            ind = i;
            d = theta - thetavals[i-1];
            dx = thetavals[i] - thetavals[i-1];
            dyv = Vel[i].Vr - Vel[i-1].Vr;
            dyt = Vel[i].Vtheta - Vel[i-1].Vtheta;
            res.Vr = Vel[i].Vr + dyv/dx*d;
            res.Vtheta = Vel[i].Vtheta + dyt/dx*d;
            break;
        }
    }
    
    return res;
}

void write_vtk(vec** bot_arr, vec** top_arr, vec** RLE, int roundpoints, int points, int* num_samples, int seednum, 
    const std::string& filename, limits limits) {
    
    std::ofstream vtk_file(filename);
    if (!vtk_file.is_open()) {
        std::cerr << "Could not open file " << filename << std::endl;
        return;
    }
    
    // Write VTK Header
    vtk_file << "# vtk DataFile Version 3.0\n";
    vtk_file << limits.length << " " << limits.planar_area << " " << limits.x_edge << " (Length, Planar Area, X limit)\n";
    vtk_file << "ASCII\n";
    vtk_file << "DATASET POLYDATA\n";
    vtk_file << "POINTS " << 2*points + roundpoints*seednum << " float\n";
   

    // Write points
    std::vector<int>* top_points =  new std::vector<int>[seednum];
    std::vector<int>* bottom_points =  new std::vector<int>[seednum];
    std::vector<int>* le_points = new std::vector<int>[seednum];
    
    int ind = 0;
    for (int i = 0; i < seednum; i++) {
        int len = num_samples[i];
        for(int j = 0; j<len; j++){
            vtk_file << top_arr[i][j].x << " " << top_arr[i][j].y << " " << top_arr[i][j].z << "\n";
            top_points[i].push_back(ind);
            ind++;
            vtk_file << bot_arr[i][j].x << " " << bot_arr[i][j].y << " " << bot_arr[i][j].z << "\n";
            bottom_points[i].push_back(ind);
            ind++;
        }
        //std::cout<<top_points[i].size()<<"\n";
        for(int k=0; k<roundpoints; k++){
            vtk_file << RLE[i][k].x << " " << RLE[i][k].y  << " " << RLE[i][k].z  << "\n";
            le_points[i].push_back(ind);
            ind++;
        }
    }
    // Do not stare into the void, for it will stare back
    std::vector<tri> triangles;
    tri temp;
    for(int row=0; row<seednum-1; row++) {
        //Back surface
        int back = top_points[row].size()-1;
        int backa = top_points[row+1].size()-1;
        int backb = top_points[row-1].size()-1;
        temp.x = top_points[row][back]; temp.y =  bottom_points[row+1][backa]; temp.z =  bottom_points[row][back]; temp.comp = 2;
        triangles.push_back(temp);
        temp.x = top_points[row][back]; temp.y = top_points[row+1][backa]; temp.z = bottom_points[row+1][backa]; temp.comp = 2;
        triangles.push_back(temp);

        ind = 0;
        if(num_samples[row] < num_samples[row+1]){
            temp.x = top_points[row][0]; temp.y = top_points[row+1][0]; temp.z = top_points[row+1][1]; temp.comp = 1;
            triangles.push_back(temp);
            temp.x = bottom_points[row][0]; temp.y = bottom_points[row+1][1]; temp.z = bottom_points[row+1][0]; temp.comp = 1;
            triangles.push_back(temp);
            ind ++;
            while(ind < num_samples[row]){
                temp.x = top_points[row][ind]; temp.y = top_points[row+1][ind]; temp.z = top_points[row+1][ind+1]; temp.comp = 1;
                triangles.push_back(temp);
                temp.x = top_points[row][ind]; temp.y = top_points[row][ind-1]; temp.z = top_points[row+1][ind]; temp.comp = 1;
                triangles.push_back(temp); 
                temp.x = bottom_points[row][ind]; temp.y = bottom_points[row+1][ind+1]; temp.z = bottom_points[row+1][ind]; temp.comp = 1;
                triangles.push_back(temp);
                temp.x = bottom_points[row][ind]; temp.y = bottom_points[row+1][ind]; temp.z = bottom_points[row][ind-1]; temp.comp = 1;
                triangles.push_back(temp);
                ind++;
            }   
        }
        else if(num_samples[row] > num_samples[row+1]){
            temp.x = top_points[row][0]; temp.y = top_points[row+1][0]; temp.z = top_points[row][1]; temp.comp = 1;
            triangles.push_back(temp);
            temp.x = bottom_points[row][0]; temp.y = bottom_points[row][1]; temp.z = bottom_points[row+1][0]; temp.comp = 1;
            triangles.push_back(temp); 
            ind+=2;
            while(ind < num_samples[row]){
                temp.x = top_points[row][ind]; temp.y = top_points[row+1][ind-2]; temp.z = top_points[row+1][ind-1]; temp.comp = 1;
                triangles.push_back(temp); 
                temp.x = top_points[row][ind]; temp.y = top_points[row][ind-1]; temp.z = top_points[row+1][ind-2]; temp.comp = 1;
                triangles.push_back(temp); 
                temp.x = bottom_points[row][ind]; temp.y = bottom_points[row+1][ind-1]; temp.z = bottom_points[row+1][ind-2]; temp.comp = 1;
                triangles.push_back(temp); 
                temp.x = bottom_points[row][ind]; temp.y = bottom_points[row+1][ind-2]; temp.z = bottom_points[row][ind-1]; temp.comp = 1;
                triangles.push_back(temp);
                ind++;
            }   
        }

        else{
            ind++;
            while(ind < num_samples[row]){
                temp.x = top_points[row][ind]; temp.y = top_points[row][ind-1]; temp.z = top_points[row+1][ind-1]; temp.comp = 1;
                triangles.push_back(temp);
                temp.x = top_points[row][ind]; temp.y = top_points[row+1][ind-1]; temp.z = top_points[row+1][ind]; temp.comp = 1;
                triangles.push_back(temp);
                temp.x = bottom_points[row][ind]; temp.y = bottom_points[row+1][ind-1]; temp.z = bottom_points[row][ind-1]; temp.comp = 1;
                triangles.push_back(temp);
                temp.x = bottom_points[row][ind]; temp.y = bottom_points[row+1][ind]; temp.z = bottom_points[row+1][ind-1]; temp.comp = 1;
                triangles.push_back(temp);
                ind++;
            }
            
        }
        
        //Leading Edge
        //if(row==50) std::cout<<triangles.size()<<std::endl;
        temp.x = top_points[row][0]; temp.y = le_points[row][0]; temp.z = le_points[row+1][0]; temp.comp = 3;
        triangles.push_back(temp); 
        temp.x = top_points[row][0]; temp.y = le_points[row+1][0]; temp.z = top_points[row+1][0]; temp.comp = 3;
        triangles.push_back(temp); 
        temp.x = le_points[row][roundpoints-1]; temp.y = bottom_points[row][0]; temp.z = bottom_points[row+1][0]; temp.comp = 3;
        triangles.push_back(temp);
        temp.x = le_points[row][roundpoints-1]; temp.y = bottom_points[row+1][0]; temp.z = le_points[row+1][roundpoints-1]; temp.comp = 3;
        triangles.push_back(temp); 
        for(int j=0; j<roundpoints-1; j++){
            temp.x = le_points[row][j]; temp.y = le_points[row][j+1]; temp.z = le_points[row+1][j+1]; temp.comp = 3;
            triangles.push_back(temp); 
            temp.x = le_points[row][j]; temp.y = le_points[row+1][j+1]; temp.z = le_points[row+1][j]; temp.comp = 3;
            triangles.push_back(temp); 
        }
    }

    //Nodules or whatever theyre called
    for(int i=0; i<roundpoints-1; i++){
        temp.comp = 2;
        temp.x = bottom_points[0][0]; temp.y = le_points[0][i+1]; temp.z = le_points[0][i];
        triangles.push_back(temp); 
        temp.x = bottom_points[seednum-1][0]; temp.y = le_points[seednum-1][i]; temp.z = le_points[seednum-1][i+1];
        triangles.push_back(temp); 
    }

    temp.x = bottom_points[0][0]; temp.y = le_points[0][0]; temp.z = top_points[0][0];
    triangles.push_back(temp);
    temp.x = bottom_points[seednum-1][0]; temp.y = top_points[seednum-1][0]; temp.z = le_points[seednum-1][0];
    triangles.push_back(temp);

    int tris = triangles.size();
    vtk_file << "POLYGONS " << tris << " " << tris*4 << "\n";
    for(int i=0; i<tris;i++) vtk_file << "3 " << triangles[i].x << " " << triangles[i].y << " " << triangles[i].z << "\n";

    vtk_file<<"CELL_DATA " << tris << "\nSCALARS Components int\nLOOKUP_TABLE default\n";
    for(int i=0; i<tris;i++) vtk_file << triangles[i].comp << "\n";

    vtk_file.close();
    std::cout << "VTK file written to " << filename << "\n";
}

int howmanypoints(std::vector<double> all[], int seednum){
    int tally = 0;
    for(int i=0; i<seednum; i++) tally+= all[i].size();
    return tally;
}

void generate3dArc(vec b1, vec b2, vec C, int num, int seednum, vec** RLE){
    vec p2 = b1; p2.y = -b1.y;
    double tol = 1e-8;
    double Rsq = b1.x*b1.x+b1.y*b1.y+b1.z*b1.z;
    double R = sqrt(R);

    //Take the dot product
    double dot = b1.x*b1.x-b1.y*b1.y+b1.z*b1.z;
    double cth = dot/Rsq;
    double th = acos(cth);

    for(int i = 0; i<num; i++){
        double numer = i+1;
        double denom = num+1;
        double thind = numer/denom*th; //Go from 1/9 to 8/9
        vec point;
        point.x = C.x + b1.x*cos(thind) + b2.x*sin(thind);
        point.y = C.y + b1.y*cos(thind) + b2.y*sin(thind);
        point.z = C.z + b1.z*cos(thind) + b2.z*sin(thind);
        RLE[seednum][i] = point;
    }
}

int main(int argc, char* argv[]){
    if(argc<2){
        std::cout<<"Usage: cwg inputs.txt (output.vtk)\n";
        return 1;
    }

    //READ INPUTS
    const double dtor = 0.01745329251994330;
    std::ifstream inpfile(argv[1]);
    std::string line = " ";
    std::string file = "";

    readinputs(inpfile, line, "beta");//Shock angle (rad)
    double beta = dtor*stod(line); file = file + line; 
    readinputs(inpfile, line, "M");//Freestream Mach number
    double M1 = stod(line); file = file + line;       
    readinputs(inpfile, line, "gamma");//Ratio of heat capacities
    double gamma = stod(line); file = file + line;    
    readinputs(inpfile, line, "h");//cone step size 
    double h = stod(line);
    readinputs(inpfile, line, "seednum");//Number of seedpoints   
    int seednum = stoi(line);   

    readinputs(inpfile, line, "span");   
    double span = stod(line);file = file + line; 
    readinputs(inpfile, line, "voffset");   
    double voffset = stod(line);file = file + line; 
    readinputs(inpfile, line, "roundpoints");   
    int roundpoints = stoi(line);

    inpfile.close();

    fs::create_directory("./Waverider Designs");
    std::string filename;
    if(argc == 2){
        filename = "./Waverider Designs/" + file + ".vtk";
    } 
    else {
        file = argv[2]; filename = "./Waverider Designs/" + file;
    }
    
    Vcone MatV = ICs(beta,M1,gamma);
    std::vector<Vcone> Vel; 
    std::vector<double> thetavals;
    double theta = beta;
    double coneangle;

    limits limits;

    while(true){
        double dV = TMS(MatV.Vr, MatV.Vtheta, gamma, theta);
        MatV.Vr += h*MatV.Vtheta;
        MatV.Vtheta += h*dV;
        theta += h;

        Vel.push_back(MatV);
        thetavals.push_back(theta);
        if(MatV.Vtheta > 0){
            coneangle = theta-h;
            break;
        }
    }
    
    //LEP Construction and Projection
    int seedsgen = (seednum+1)/2; //Number of seed points to generate (we only do half)
    vec* LEP_points = new vec[seednum];
    double xpos = -span/2;
    double rpl = 0.064/2; //Payload Radius
    double roffset = 1.5;
    double dip;

    //Option 1
    //dip = rpl*roffset-rpl*exp(-0.5*xpos*xpos/rpl/rpl);
    //z_center = coneslope*rpl*(roffset-1);
    auto LEPslope = [](double xpos, double rpl){return xpos/rpl*exp(-xpos*xpos/2/rpl/rpl);};

    //Option 2
    //dip = xpos*xpos/2/rpl+10*rpl;
    //z_center = coneslope*10*rpl;
    //auto LEPslope = [](double xpos, double rpl){return xpos/rpl;};

    double coneslope = 1/tan(beta);
    double conesurface;
    for(int i=0; i<seednum; i++){
        dip = rpl*roffset-rpl*exp(-0.5*xpos*xpos/rpl/rpl);
        LEP_points[i].x = xpos;
        LEP_points[i].z = coneslope*sqrt(LEP_points[i].x*LEP_points[i].x+dip*dip);
        LEP_points[i].y = sqrt(LEP_points[i].z*LEP_points[i].z*tan(beta)*tan(beta) - xpos*xpos);
        xpos+=span/(seednum-1);
    }
    double z_center = coneslope*rpl*(roffset-1);
    double y_center = z_center*tan(beta);

    double zmax = LEP_points[0].z;

    //Surface Streamline Interpolation and Tracing
    double SSF = 2; //Step safety factor (Dubbed by the honorable Kenshiro Lim)
    int mesh_res = ceil(seednum / 2);
    double spline_step = (zmax-z_center) / (SSF*mesh_res);

    // number of samples corresponding to each seed point
    int* samplenum = new int[seednum]; 
    double* all_rho = new double[seednum];
    std::vector<double>* all_theta = new std::vector<double>[seednum];
    std::vector<double>* all_r = new std::vector<double>[seednum];

    int num_points = 0;
    int points_in_row;
    for(int i=0; i<seednum; i++){
        all_rho[i] = atan2(LEP_points[i].y,LEP_points[i].x);
        if(i<seednum/2) points_in_row = i+1;
        else points_in_row = seednum-i;
        samplenum[i] = points_in_row;
        num_points += points_in_row;
    }

    //fill in theta and r arrays
    for(int i=0; i<seednum; i++){
        double ith_r = sqrt(LEP_points[i].x*LEP_points[i].x+LEP_points[i].y*LEP_points[i].y+LEP_points[i].z*LEP_points[i].z);
        double ith_theta = atan(sqrt(LEP_points[i].x*LEP_points[i].x + LEP_points[i].y*LEP_points[i].y) / LEP_points[i].z);

        double curr_r = ith_r;
        double curr_theta = ith_theta;
        double curr_z = curr_r * cos(curr_theta);

        double prev_r;
        double prev_theta;
        double prev_z;
        double final_r;

        std::vector<double> spline_r;
        std::vector<double> spline_theta;
        spline_r.push_back(curr_r);
        spline_theta.push_back(curr_theta);
        
        while(true){
            Vcone V_curr = Interp(curr_theta, Vel, thetavals);
            prev_r = curr_r;
            prev_theta = curr_theta;
            prev_z = prev_r * cos(prev_theta);

            // Step the radial and theta components
            curr_r = curr_r + V_curr.Vr * spline_step;
            curr_theta = curr_theta + V_curr.Vtheta * spline_step / curr_r; 
            curr_z = curr_r * cos(curr_theta);
      
            spline_r.push_back(curr_r);
            spline_theta.push_back(curr_theta);

            if(curr_z >= zmax){
                //Interpolation factor
                double factor = (zmax - prev_z)/(curr_z-prev_z);
                final_r = prev_r + factor*(curr_r-prev_r);
                break;
            }
        }
        
        int num_points = samplenum[i]; 
        all_r[i].push_back(ith_r);
        all_theta[i].push_back(ith_theta);
        
        if(num_points>1){
            //int sample_ind = 0;
            for(int j=1; j<num_points; j++){
                double r_sample  = (j)*(final_r-ith_r)/(num_points-1)+ith_r;
                int r_ind = 0;
                while(spline_r[r_ind]<=r_sample && r_ind < spline_r.size()) r_ind++;
                //if(r_ind == 0) std::cout<<ith_r<<" "<<r_sample<<" "<<spline_r[0]<<"\n";

                double theta_inter = (r_sample - spline_r[r_ind-1])/(spline_r[r_ind]-spline_r[r_ind-1]);
                double theta_sample = spline_theta[r_ind-1] + theta_inter*(spline_theta[r_ind]-spline_theta[r_ind-1]);

                all_r[i].push_back(r_sample);
                all_theta[i].push_back(theta_sample);
                //sample_ind++;
            }
        }
        //std::cout<<all_r[i].size()<<"\n";
    }

    //Cylindrical to Cartesian
    vec** top_points = new vec*[seednum];
    vec** bot_points = new vec*[seednum];
    for(int i=0; i<seednum; i++) {
        top_points[i] = new vec[samplenum[i]];
        bot_points[i] = new vec[samplenum[i]];
    }
   
    //Bottom Surface
    for(int i=0; i<seednum; i++){
        for(int j=0; j<samplenum[i]; j++){
            bot_points[i][j].x = all_r[i][j]*sin(all_theta[i][j])*cos(all_rho[i]);
            bot_points[i][j].y = all_r[i][j]*sin(all_theta[i][j])*sin(all_rho[i])-y_center;
            bot_points[i][j].z = all_r[i][j]*cos(all_theta[i][j])-z_center;
        }
    }
    
    for(int i=seedsgen; i<seednum; i++){
        int mirror = seednum-1-i;
        for(int j=0; j<samplenum[mirror]; j++){
            bot_points[i][j].x = -bot_points[mirror][j].x;
            bot_points[i][j].y = bot_points[mirror][j].y;
            bot_points[i][j].z = bot_points[mirror][j].z;
        }
    }
    
    //Top Surface
    for(int i=0; i<seednum; i++){
        int len = samplenum[i];
        double xpos = LEP_points[i].x;
        double ypos = LEP_points[i].y-voffset;
        double zstart = LEP_points[i].z;
        for(int j=0; j<len; j++){
            top_points[i][j].x = xpos;
            top_points[i][j].y = ypos-y_center;
            if(len == 1) top_points[i][j].z = zstart-z_center;
            else top_points[i][j].z = j*(zmax-zstart)/(len-1)+zstart-z_center;
        }
    }

    //Rounded leading edge hack
    vec** roundededge = new vec*[seednum];
    for(int i = 0; i<seednum; i++) roundededge[i] = new vec[roundpoints];
    double tol = 1e-5; //Tolerance for how close to (0,1) our inward vector is

    //Endcases
    //Slope of back
    vec temp; temp.z=0; //IGotta stop calling these variables temp
    double vx = top_points[0][0].x-top_points[1][1].x; //we know this actually
    double vy = top_points[0][0].y-top_points[1][1].y;
    vec norm;
    double factor = voffset/2/vx; //remember that voffset gets subtracted
    vx*=factor; vy*=factor; norm.x = -vy; norm.y = vx; norm.z = 0;
    temp.x = -vx; temp.y = -vy;
    vec center;
    center.x = top_points[0][0].x + norm.x;
    center.y = top_points[0][0].y + norm.y;
    center.z = top_points[0][0].z + norm.z;
    vec b1; b1.x = vy; b1.y = -vx; b1.z = 0;
    generate3dArc(b1,temp,center,roundpoints,0,roundededge);
    limits.x_edge = sqrt(vy*vy+vx*vx)+fabs(center.x);

    for(int j=0; j<roundpoints; j++) {
            roundededge[seednum-1][j].x = -roundededge[0][j].x;
            roundededge[seednum-1][j].y = roundededge[0][j].y;
            roundededge[seednum-1][j].z = roundededge[0][j].z;
    }

    //Black magic (DO NOT LET ME COOK)
    for(int i = 1; i<seednum-1; i++){
        // No longer normal to leading edge, simply normal to the shock cone
        double xpos = top_points[i][0].x;
        double slope = LEPslope(xpos,rpl);
        double dx = -slope; //This is the component of the vector (dx,1) in the xz-plane that points inward normal to the LEP
        
        //Use plane containing point directly behind and point directly behind and LEP point directly outside
        double bruh = fabs(dx);
        if(fabs(top_points[i+1][1].y - top_points[i][0].y)<tol) {
            temp.x = 0; temp.y = 0; temp.z = 1;
        }
        else{
            //look to the streamline one point inward.
            double xdiff,ydiff;
            xdiff = top_points[i+1][1].x - top_points[i][0].x;
            ydiff = top_points[i+1][1].y - top_points[i][0].y;
            //figure how much we need to scale our vector (dx,1) to get to the next streamline
            double scale = fabs(xdiff/dx);

            //So our vector tangent along the surface is (xdiff, ydiff, scale)
            temp.x = xdiff; temp.y = ydiff; temp.z = scale;
        }
        //We need to rotate this vector 90 deg down (if this works, convert to triple cross product formula)
        vec down; down.x = 0; down.y = 1; down.z = 0; // y is positive because of how the waverider is generated

        //inward = temp x (down x temp)
        //a x (b x a) = (a.a)b - (a.b)a
        double aa = temp.x*temp.x+temp.y*temp.y+temp.z*temp.z;
        double ab = -temp.y;
        vec inward; inward.x = ab*temp.x; inward.y = aa+ab*temp.y; inward.z = ab*temp.z;

        //Find how much we need to scale our "inward" vector now
        double ygoal = voffset/2;
        factor = ygoal/inward.y;
        
        vec b1; b1.x = -factor*inward.x; b1.y = -ygoal; b1.z = -factor*inward.z; //This is our first basis vector
        double magtemp = sqrt(temp.x*temp.x+temp.y*temp.y+temp.z*temp.z);
        double magb1 = sqrt(b1.x*b1.x+b1.y*b1.y+b1.z*b1.z);
        temp.x *= -magb1/magtemp; temp.y *= -magb1/magtemp; temp.z *= -magb1/magtemp;
        
        
        vec b2 = temp;
        vec center; center.x = top_points[i][0].x-b1.x; center.y = top_points[i][0].y-b1.y; center.z = top_points[i][0].z-b1.z;
        generate3dArc(b1,temp,center,roundpoints,i,roundededge);

        //std::cout<<v1.x<<" "<<v1.y<<" "<<v1.z<<" / "<<v2.x<<" "<<v2.y<<" "<<v2.z<<" / "<<nvec.x<<" "<<nvec.y<<" "<<nvec.z<<" \n ";
        //std::cout<< slope << "\n";
        if(i>=seednum/2) {//Just put the fries in the bag
            for(int j=0; j<roundpoints; j++) {
                roundededge[i][j].x = -roundededge[seednum-1-i][j].x;
                roundededge[i][j].y = roundededge[seednum-1-i][j].y;
                roundededge[i][j].z = roundededge[seednum-1-i][j].z;
            }
        }
    }

    write_vtk(bot_points, top_points, roundededge, roundpoints, num_points, samplenum, seednum, filename, limits);
    return 0;
}