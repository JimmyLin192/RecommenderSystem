function [predict_label, accuracy] = oneClassSVM ()

filename = './sampleData/joint.csv';
[labels, instances] = libsvmread(filename)

nInsts = size(labels);
nTrainInst = 100;

train_label = labels(1:nTrainInst, :);
test_label = labels(nTrainInst:nInsts, :);
train_instance = instances(1:nTrainInst, :);
test_instance = instances(nTrainInst:nInsts,:);

model = svmtrain(train_label, train_instance, '-s 2');

[predict_label, accuracy, prob_estimates] = svmpredict(test_label, test_instance, model);

end