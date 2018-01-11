\name{drs}
\alias{drs}
\docType{data}
\title{Diabetic Retinopathy Study (DRS)}
\description{
The Diabetic Retinopathy Study (DRS) was performed to determine whether the
onset of blindness in 197 high-risk diabetic patients could be delayed by
laser treatment. The treatment was administered to one
randomly-selected eye in each patient, leaving the other eye untreated. Thus,
there are 394 observations, which are clustered by patient since the level of
risk will tend to vary between patients. A failure occurred when visual
acuity dropped to below 5/200. All patients had a visual acuity of at least
20/100 at the beginning of the study.
}
\usage{
data("drs")
}
\format{A data frame with 394 rows and 8 columns. There are two rows for each
subject, one row for each eye:
\describe{
  \item{subject_id}{unique identifier for each subject}
  \item{eye}{subject's eye, where 1=right and 2=left}
  \item{time}{the observed follow-up time}
  \item{status}{outcome at the end of the observation period, where
  1=blindness and 0 indicates censorship}
  \item{treated}{a binary covariate, where 1=treated or 0=untreated}
  \item{age_at_onset}{age (in years) at the onset of diabetes}
  \item{laser_type}{type of laser used for treatment, where 1=xenon, 2=argon}
  \item{diabetes_type}{type of diabetes, where 1=juvenile (age at dx < 20) and 2=adult}
}}
\source{
\url{http://www.mayo.edu/research/documents/diabeteshtml/doc-10027460}
}
\examples{
\dontrun{
data(drs)

# Clustered by subject
fit.drs <- fitfrail(Surv(time, status) ~ treated + cluster(subject_id), 
                    drs, frailty="gamma")

fit.drs

# Variance estimates
vcov(fit.drs)

# Plot the estimated cumulative baseline hazard
plot(fit.drs, type="cumhaz")
}
}
\keyword{datasets}

