
% Strengths of the exchange field and spin-orbit field
h = [0 1/10 1/5 1 5 10];
a = [0 1/10 1/5 1 5 10];

% Initiate the pool of parallel workers
matlabpool local 12

iterations = [length(h), length(a)];
parfor n=1:prod(iterations)
    % Retrieve the indices corresponding to this iteration and
    [i,j] = ind2sub(iterations,n);
    
    % Information about parallel execution
    task   = getCurrentTask;
    taskID = task.ID;
    
    % Display progress information
    disp(sprintf('-- Worker %2.0f: commencing calculation for h=%2.2f, a=%2.2f', taskID, h(i), a(j)));
    
    % Instantiate a ferromagnet
    f = Ferromagnet( linspace(0,1,128), linspace(0.01, 1.3 , 20) );

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
    filename = sprintf('output/ferromagnet_h%2.2f_a%2.2f.mat', h(i), a(j))
    save(filename, 'f');
    
    % Display progress information
    disp(sprintf('-- Worker %2.0f: saving results to "%s".', taskID, filename));
end