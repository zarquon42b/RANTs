% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/landmarkTransform.R
\name{transform2antsTrafo}
\alias{transform2antsTrafo}
\title{convert an affine transform from matrix to antsrTransform}
\usage{
transform2antsTrafo(affinemat, IJK2RAS = diag(c(-1, -1, 1,
  1)[seq_len(nrow(affinemat))]))
}
\arguments{
\item{affinemat}{a 4x4 (3D case) or 3x3 (2D case) homogenous transform matrix}

\item{IJK2RAS}{transform from point to image space.}

\item{file}{filename}
}
\value{
returns object of class antsrTransform
}
\description{
convert an affine transform from matrix to antsrTransform (implementation of itk::Transform
}
