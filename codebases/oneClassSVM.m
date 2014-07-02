function [predict_label, accuracy] = oneClassSVM ()

train_filename = './sampleData/joint_train.csv';
test_filename = './sampleData/joint_test.csv';
[train_labels, train_instances] = libsvmread(train_filename);
[test_labels, test_instances] = libsvmread(test_filename);


% one-class support vector machine
model = svmtrain(train_labels, train_instances, '-s 2 -t 4');

[predict_label, accuracy, prob_estimates] = svmpredict(test_labels, test_instances, model);

predict_label
accuracy

end
