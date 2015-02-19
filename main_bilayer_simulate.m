function simulate_bilayer_bulk()
    % Log output to a file
    mkdir output
    diary output/simulation.log
    
    % How many positions and energies to consider
    positions = 32;
    energies  = 32;

    % Strengths of the exchange field, strength of spin-orbit field, length
    h = [0 1/10  1 10 100 1000];
    a = [0 1/100 1/10 1 10 100];
    d = [1/10000 1/1000  1/100 1/10 1 10 100 1000 10000];


    iterations = [length(d), length(h), length(a)];
    disp('Parallel execution enabled.');
    parfor n=1:prod(iterations)
%      disp(':: Parallel execution disabled.');
%      for n=1:prod(iterations)
        % Retrieve the indices corresponding to this iteration and
        [k,i,j] = ind2sub(iterations,n);

        % Information about parallel execution
        task = getCurrentTask();
        if isempty(task)
            taskID = 0;
        else
            taskID = task.ID;
        end

        % Display progress information and make sure results can be saved
        fprintf(':: Worker %2.0f: commencing calculation for d=%2.4f, h=%2.4f, a=%2.4f.\n', taskID, d(k), h(i), a(j));
        filename = sprintf('ferromagnet_d%2.4f_h%2.4f_a%2.4f', d(k), h(i), a(j));
        parsave(filename, 'f');
        pause(0.25);

        % Instantiate a ferromagnet
        f = Ferromagnet( d(k)*linspace(0,1,positions), linspace(0.05, 1.5, energies) );

        % Change the exchange and spin-orbit fields of the ferromagnet
        f.exchange  = [h(i),0,0];
        f.spinorbit = SpinVector.RashbaDresselhaus(a(j),pi/4);

        % Set the boundary conditions for the ferromagnet (i.e. connect it to a bulk superconductor)
        f.interface_left = 1;
        for m=1:length(f.energies)
            f.boundary_left(m) = Superconductor.Bulk(f.energies(m), 1.0);
        end

        % Gradually increase the precision
        for N=([4 8 32]*length(f.positions))
            % Update the maximum grid size and error tolerance
            f.sim_grid_size = N;
            f.sim_error_abs = 0.05;
            f.sim_error_rel = 0.05;
            
            % Update the internal state of the ferromagnet
            f.update;
            
            % Store the results of the simulation 
            parsave(filename, f);

            % Display progress information
            fprintf(':: Worker %2.0f: saving results to "%s" (grid size:%3.f)\n', taskID, filename, N);
            pause(1);
        end
    end
    
    % Turn off the logging function
    diary off
end

function parsave(filename, f)
    save(sprintf('output/%s.mat', filename), 'f');
end