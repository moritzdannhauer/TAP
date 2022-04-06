function nvec=trinormal(tri,p)
% Function to compute normal vector to triangle
% Input:    tri <triangular index matrix>
%           p nx3 array of vertices
% Output:   nvec <nx3> array of normal vector
% JOK180609
% Version: 1

% I/O check
% to be completed
% Check polarity of elements

% Construct vectors 
v = [p(tri(:,3),1)-p(tri(:,1),1), p(tri(:,3),2)-p(tri(:,1),2), p(tri(:,3),3)-p(tri(:,1),3)];
w = [p(tri(:,2),1)-p(tri(:,1),1), p(tri(:,2),2)-p(tri(:,1),2), p(tri(:,2),3)-p(tri(:,1),3)];
% Calculate cross product
normvec = [v(:,2).*w(:,3)-v(:,3).*w(:,2), ...
    -(v(:,1).*w(:,3)-v(:,3).*w(:,1)), ...
    v(:,1).*w(:,2)-v(:,2).*w(:,1)];
% Normalize
lnvec = sqrt(sum(normvec.*normvec,2));
nvec = normvec./repmat(lnvec,1,3);