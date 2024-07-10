function [source_rot, cav, curlz] = computeRotationalSouces(p,dx,dy,evaluation_points,X,Y)

    cav = nan(size(p)); curlz = cav;
    
    % curl calculation
    for ii = 1:numel(p(1,1,:))
        dx(:,:,ii) = inpaint_nans(dx(:,:,ii));
        dy(:,:,ii) = inpaint_nans(dy(:,:,ii));
        
        %cav(:,:,ii) = curl(X,Y,dx(:,:,evaluation_points(ii)), dy(:,:,evaluation_points(ii)));
        [curlz(:,:,ii), cav(:,:,ii)] = curl(X,Y,dx(:,:,ii), dy(:,:,ii));
        % cav(:,:,ii) = inpaint_nans(cav(:,:,ii)); curlz(:,:,ii) = inpaint_nans(curlz(:,:,ii));
    end

    % source point
    source_rot = nan( 2, length(evaluation_points) );
    for ii = 1:length(evaluation_points)    
        [yy,xx] = find( curlz(:,:,ii) == max( reshape(curlz(:,:,ii), 1, [] ) ) );
        if numel(yy) == 1
            source_rot(1,ii) = xx; source_rot(2,ii) = yy;
        else
            source_rot(1,ii) = NaN; source_rot(2,ii) = NaN; 
        end
    end

end