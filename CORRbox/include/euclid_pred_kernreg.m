function prediction = euclid_pred_kernreg(C3dnew, Cloc, theModel, nameFunction)

% checker
if (~corraux_checker(C3dnew))
    errmsg = "* ";
    errmsg = errmsg + nameFunction + " : an input is not a proper object. Use 'corr_init'.";
    error(errmsg);
end

C3data = C3dnew.data;
if (size(C3data,1)~=size(Cloc,1))
    errmsg = "* ";
    errmsg = errmsg + nameFunction + " : the provided data does not match with the trained model dimension.";
    error(errmsg);
end

% find the localized coordinates for the new data
X_test = euclid_localise(Cloc, C3data);

% predict
prediction = predict(theModel, X_test);


end