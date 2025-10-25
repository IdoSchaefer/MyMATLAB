function result = invHessRosen(x)
    result = [0.5,   x(1);
              x(1),  3*x(1)^2 - x(2) + 0.005]/(200*x(1)^2 - 200*x(2) + 1);
end