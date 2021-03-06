% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/landmarkTransform.R
\name{transform2mat}
\alias{transform2mat}
\title{write an affine transform matrix to ITK mat format}
\usage{
transform2mat(affinemat, file, IJK2RAS = diag(c(-1, -1, 1,
  1)[seq_len(nrow(affinemat))]))
}
\arguments{
\item{affinemat}{a 4x4 (3D case) or 3x3 (2D case) homogenous transform matrix}

\item{file}{filename}

\item{IJK2RAS}{transform from point to image space.}
}
\value{
returns the character of the filename
}
\description{
write an affine transform matrix to ITK mat format
}
