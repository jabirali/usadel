
% Strengths of the exchange field and spin-orbit field
h = [0 1/10 1/5 1 5 10];
a = [0 1/10 1/5 1 5 10];

% Create a storage place for the results of simulations
singlets = State.empty(length(h), length(a), length(f.positions), length(f.energies));

% Initiate the pool of parallel workers
%matlabpool local 12

iterations = [length(h), length(a)];
%parfor
for n=1:prod(iterations)
    % Retrieve the indices corresponding to this iteration
    [i,j] = ind2sub(iterations,n);
    
    % Instantiate a ferromagnet
    f = Ferromagnet( linspace(0,1,5), 1:0.2:1.3 );

    % Change the exchange and spin-orbit fields of the ferromagnet
    f.exchange  = [h(i),0,0];
    f.spinorbit = SpinVector.RashbaDresselhaus(a(j),pi/4);

    % Set the boundary conditions for the ferromagnet (i.e. connect it to a bulk superconductor)
    f.interface_left = 1;
    for m=1:length(f.energies)
        f.boundary_left(m) = Superconductor.Bulk(f.energies(m), 1.0);
    end
      
    % Update the internal state of the ferromagnet
    f.update;
        
    % Store the results of the simulation 
    filename = ['ferromagnet_h' num2str(h(i)) '_a' num2str(a(j)) '.mat']
    save filename f;
end