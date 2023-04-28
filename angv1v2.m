function angle = angv1v2(v1,v2)

angle = acosd( dot(v1,v2) / (norm(v1)*norm(v2)) );

end