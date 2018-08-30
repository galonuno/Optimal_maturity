function c_out=c_path(c_in,rss,steady,paths,parameters)
path_out=compute_path(c_in,rss,steady,paths,parameters);
c_out=path_out.c_out;
% disp(iter)
% disp(residual)