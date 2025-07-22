function [ p ] = recon_dat_store( p )
% For loading positions from a previously stored recostruction

for ii = 1:p.numscans

    switch p.scan.type
        case 'custom'
            pos_file = p.scan.custom_positions_source; % added by Arthur B
            if exist(pos_file,'file')
                psn = load(pos_file);
                %{
                Note that in ptycho\+core\initialize_ptycho.m
                	line 94: p = scans.read_positions(p) -> ulitmately leads to
                this file:
                	line 103: p = core.ptycho_adjust_positions(p);
                calls ptycho\+core\ptycho_adjust_positions.m
                which on the first line execs:
                	p.positions_real = -p.positions_real;
                then performs:
                	% Convert to pixels
                	p.positions = p.positions_real./p.dx_spec;
                
                so we need to multiply by p.dx_spec and multiply by -1

                %}
                p.positions_real = -1*psn.recon_dat.shelves_op.positions* ...
                                    p.dx_spec; % 
            else
                error('Could not find function or data file %s', pos_file);
            end
            
        otherwise
            error('Unknown scan type %s.', p.scan.type);
    end
    
    p.numpts(ii) = size(positions_real,1);
    p.positions_real = [p.positions_real ; positions_real]; %append position - 

end


end

