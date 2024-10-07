classdef getSparsityLib
    methods(Static)
        function S = mass_transport_1D(Nx,Ny)
            totalSize = Nx;
            e = ones(totalSize,1);
            A = spdiags([e e e],-1:1,totalSize,totalSize);
            S = full(A);
        end
        function S = cahn_hilliard_2D(Nx,Ny)
            e = ones(Ny,1);
            A = spdiags([e e e e e],-2:2,Ny,Ny);
            v = full(A);
            S = repmat(v,Nx,Nx);
        end
        function S = cahn_hilliard_1D(Nx,Ny)
            %totalSize = Nx;
            e = ones(Nx,1);
            A = spdiags([e e e e e],-2:2,Nx,Nx);
            v = full(A);
            %S = ones(totalSize,totalSize);
            S = v;
        end
    end
end
