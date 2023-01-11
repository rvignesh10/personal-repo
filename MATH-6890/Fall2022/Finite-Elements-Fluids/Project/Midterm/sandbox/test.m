%%
files = dir("../solveFiles/*.txt");

for i=1:length(files)
    str = strcat(files(i).folder,'/');
    str = strcat(str, files(i).name);
    fname = extractBefore(str,'.txt');
    fname = extractAfter(fname,'Files/');
    m = readmatrix(str);
    FE_Data.(fname) = m;
end

nodes = FE_Data.FiniteElementSpace(3,2);
S = sparse(FE_Data.Stabilization_0(:,1), FE_Data.Stabilization_0(:,2), FE_Data.Stabilization_0(:,3), nodes, nodes);