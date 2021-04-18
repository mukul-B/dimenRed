#ifndef PATERNREC_DISCRIMINANTCASES_H
#define PATERNREC_DISCRIMINANTCASES_H

#endif //PATERNREC_DISCRIMINANTCASES_H

// base class ClassifierClass is created , inherited by each classifier case class
// it hold the values of matrix mu , covariance matrix and priors
// methods:
// virtual getDecision method is defined , which help in virtualization of this function
//getErrorBound -  given input of Beta ,it returns the theoretic error bond
//getEucldianDistance - for given input sample , it return the difference from the mu's of two sample

class ClassifierClass {

    Matrix m1, cv1;
    Matrix m2, cv2;
    double prior1, prior2;

public:
    ClassifierClass() {}

    ClassifierClass(Matrix mu1, Matrix mu2,
                    Matrix cov1, Matrix cov2, double prior1, double prior2) {
        this->m1 = mu1;
        this->m2 = mu2;
        this->cv1 = cov1;
        this->cv2 = cov2;
        this->prior1 = prior1;
        this->prior2 = prior2;
    }

    virtual double getDecision(Matrix x) {
        return 0;
    };

    double getEucldianDistance(Matrix x) {

        double d1 = (-1) * (((x - m1).trans()) * (x - m1)).getDouble();
        double d2 = (-1) * (((x - m2).trans()) * (x - m2)).getDouble();
        return d1 - d2;
    };

    double getErrorBound(double beta) {
        double kOfBeta = ((((m1 - m2).trans()) * (beta * (1 - beta) * (0.5))) *
                          (((cv1 * (1 - beta)) + (cv2 * beta)).inv() * (m1 - m2))).getDouble()
                         + (0.5) * log(((cv1 * (1 - beta)) + (cv2 * beta)).deter() /
                                       (pow(cv1.deter(), 1 - beta) * pow(cv2.deter(), (beta))));
        double errorBound = pow(prior1, beta) * pow(prior2, (1 - beta)) * pow(exp(1), (-1) * kOfBeta);
        return errorBound;
    }

};

//EuclideanClass use the getEucldianDistance of base class to return decision value
class EuclideanClass
        : public ClassifierClass {
public:
    EuclideanClass(Matrix m1, Matrix m2, Matrix cov1, Matrix cov2, double prior1, double prior2)
            : ClassifierClass(m1, m2, cov1, cov2, prior1, prior2) {}

    double getDecision(Matrix x) {
        return getEucldianDistance(x);
    }
};

// Case1 Classifier , inheriting ClassifierClass
// it calculate x0 and wt for the decision boundary
// override getDecision of base class , using x0 and wt

class Case1
        : public ClassifierClass {
    Matrix x0, wt;
public:
    Case1() {}

    Case1(Matrix m1, Matrix m2, Matrix cov1, Matrix cov2, double prior1, double prior2)
            : ClassifierClass(m1, m2, cov1, cov2, prior1, prior2) {
        double sigma = 0;
        if (cov1 == cov2) {
            sigma = pow(cov2.getSigmaSqr(), 0.5);
        }
        if (sigma == 0) {
            cout << "Case1 does'nt apply here as coviarance condition is not satisfied'\n";

        } else
            cout << "Case1 classifier with sigma as : " << sigma << "\n";

        x0 = ((m1 + m2)) * (0.5);
        wt = (m1 - m2).trans();
        double cof = (sigma * sigma) / (((m1 - m2).trans()) * (m1 - m2)).getDouble();
        double probilityRatio = log(prior1 / prior2);
        x0 = x0 - ((m1 - m2) * (cof * probilityRatio));
    }

    double getDecision(Matrix x) {

        Matrix decision = (wt * (x - x0));
        double discriminate_value = decision.getDouble();
        return discriminate_value;
    }
};

// Case3 Classifier , inheriting ClassifierClass
// it calculate Wi ,wi and w0i for the classes
// override getDecision of base class , deriving decision values from case3 equation
class Case3
        : public ClassifierClass {

    Matrix W1, w1;
    double w01;
    Matrix W2, w2;
    double w02;
public:
    Case3() {}

    Case3(Matrix m1, Matrix m2, Matrix cov1, Matrix cov2, double prior1, double prior2)
            : ClassifierClass(m1, m2, cov1, cov2, prior1, prior2) {

        W1 = (cov1.inv()) * (-0.5);
        w1 = (cov1.inv()) * m1;
        w01 = ((m1.trans() * (cov1.inv() * m1)) * (-0.5)).getDouble()
              - (log(cov1.deter()) * (0.5))
              + log(prior1);

        W2 = (cov2.inv()) * (-0.5);
        w2 = (cov2.inv()) * m2;
        w02 = ((m2.trans() * (cov2.inv() * m2)) * (-0.5)).getDouble()
              - (log(cov2.deter()) * (0.5))
              + log(prior2);
    }

    double getDecision(Matrix x) {

        double decision1 = ((x.trans()) * (W1 * x) + (w1.trans()) * x).getDouble() + w01;
        double decision2 = ((x.trans()) * (W2 * x) + (w2.trans()) * x).getDouble() + w02;
        double discriminate_value = decision1 - decision2;
        return discriminate_value;
    }
};

// treshold based classifier

class ThresholdBased {

    Matrix W1;
    Matrix mu;
    double w01;

public:
    ThresholdBased() {}

    ThresholdBased(Matrix m1, Matrix cov1) {
        W1 = (cov1.inv()) * (-0.5);
        mu = m1;
        w01 = 1 / (pow(cov1.deter(), 0.5) * atan(1) * 4 * 2);
    }

    double getDecision(Matrix x) {
        double decision1 =  exp((((x - mu).trans()) * (W1 * (x - mu))).getDouble());
        return decision1;
    }
};
