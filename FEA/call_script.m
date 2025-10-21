clear;
vtkFile = 'vsurf00010000.vtk';
inputFile = 'input (1).sdf';
[P_inf, T_inf, M_inf, rho_inf, pran, y, Rgas] = load_input(inputFile);
[T_e, P_e, U_e, V_e, W_e, centroids, norms, areas] = load_vtk_surf(vtkFile);

