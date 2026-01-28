function Hankel = H (L,z)

  N = size (z,1);
  m = size (z,2);
  Hankel = zeros (L*m,N-L+1);   % Initialize with zeros of the correct size

  for i = 1:N-L+1
    for j = 0:L-1
      Hankel(j*m+1:j*m+m,i) = z(i+j,:)';
    end
  end

end

