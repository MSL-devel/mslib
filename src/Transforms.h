/*
----------------------------------------------------------------------------
This file is part of MSL (Molecular Software Libraries)
 Copyright (C) 2010 Dan Kulp, Alessandro Senes, Jason Donald, Brett Hannigan,
 Sabareesh Subramaniam, Ben Mueller

This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, 
 USA, or go to http://www.gnu.org/copyleft/lesser.txt.
----------------------------------------------------------------------------
*/

#ifndef TRANSFORMS_H
#define TRANSFORMS_H

#include <set>

#include "AtomPointerVector.h"
#include "Residue.h"
#include "Quaternion.h"
#include "SphericalPoint.h"
#include <math.h>


namespace MSL { 
class Transforms {

	public:
		Transforms();
		Transforms(const Transforms & _transform);
		~Transforms();

		/***************************************************************************************
		 *
		 *  TRANSFORMATIONS APPLIED TO A SINGLE ATOM
		 *
		 ***************************************************************************************/
		// translation
		void translate(Atom & _atom, const CartesianPoint & _p);

		// rotation
		void Xrotate(Atom & _atom, double _degrees);
		void Yrotate(Atom & _atom, double _degrees);
		void Zrotate(Atom & _atom, double _degrees);
		void rotate(Atom & _atom, double _degrees, const CartesianPoint & _axisFromRotCenter, const CartesianPoint & _rotCenter=CartesianPoint(0.0, 0.0, 0.0));
		void rotate(Atom & _atom, const Matrix & _rotMatrix, const CartesianPoint & _rotCenter=CartesianPoint(0.0, 0.0, 0.0));

		/*****************************************************
		 *  Align rotates the atom so that it points in the
		 *  same direction of the target (i.e. the angle
		 *  _atom - _rotCenter - _target is = 0)
		 *
		 *  Orient operates a torsional rotation to orient the
		 *  _atom in the same direction of the _target with
		 *  respect to the axis of rotation defined bv _axis1 - _axis2 
		 *  (i.e. the dihedral(_atom, _axis1, _axis2, _target = 0)
		 *****************************************************/
		bool align(Atom & _atom, const CartesianPoint & _target, const CartesianPoint & _rotCenter=CartesianPoint(0.0, 0.0, 0.0));
		bool orient(Atom & _atom, const CartesianPoint & _target, const CartesianPoint & _axis1, const CartesianPoint & _axis2);

		// Internal functions of above, but also useful to be public.
		bool align(CartesianPoint & _object, const CartesianPoint & _target, const CartesianPoint & _rotCenter=CartesianPoint(0.0, 0.0, 0.0));
		bool orient(CartesianPoint & _object, const CartesianPoint & _target, const CartesianPoint & _axis1, const CartesianPoint & _axis2);


		/***************************************************************************************
		 *
		 *  TRANSFORMATIONS APPLIED TO AN ATOM VECTOR
		 *
		 ***************************************************************************************/
		void translate(AtomPointerVector & _atoms, CartesianPoint _p);
		void Xrotate(AtomPointerVector & _atoms, double _degrees);
		void Yrotate(AtomPointerVector & _atoms, double _degrees);
		void Zrotate(AtomPointerVector & _atoms, double _degrees);
		void rotate(AtomPointerVector & _atoms, double _degrees, const CartesianPoint & _axisFromRotCenter, const CartesianPoint & _rotCenter=CartesianPoint(0.0, 0.0, 0.0));
		void rotate(AtomPointerVector & _atoms, const Matrix & _rotMatrix, const CartesianPoint & _rotCenter=CartesianPoint(0.0, 0.0, 0.0));
		
		/*****************************************************
		 *  For the atom std::vector, the align and orient operations 
		 *  are applied to an external point _reference, and the
		 *  _atoms are moved according to the same transformation
		 *****************************************************/
		bool align(AtomPointerVector & _atoms, const CartesianPoint & _reference, const CartesianPoint & _target, const CartesianPoint & _rotCenter=CartesianPoint(0.0, 0.0, 0.0));
		bool orient(AtomPointerVector & _atoms, const CartesianPoint & _reference, const CartesianPoint & _target, const CartesianPoint & _axis1, const CartesianPoint & _axis2);


		/*******************************************************
		 *  For two atom vectors, align atoms via moving geometric centers on top, 
                 *  solving for a optimal rotation matrix.  Return value is a bool
                 *  which is true when alignment is proper.
		 *******************************************************/
		bool rmsdAlignment(AtomPointerVector &_align, AtomPointerVector &_ref);
		bool rmsdAlignment(AtomPointerVector &_align, AtomPointerVector &_ref, AtomPointerVector &_moveable);
		/*
		double align(std::vector<Residue *> &_align, std::vector<Residue *> &_ref); // Align backbone atoms of each residue
		double align(std::vector<Residue *> &_align, std::vector<Residue *> &_ref,std::vector<Residue *> &_moveable); // Align backbone atoms of each residue
		double align(std::vector<Residue *> &_align, std::vector<Residue *> &_ref,AtomPointerVector &_moveable); // Align backbone atoms of each residue
		*/

		/*******************************************************
		 * DIRECT DISTANCE AND ANGLE EDITING IN A SYSTEM
		 * if the bonding information is stored in the atoms
		 * the last atom and all those connected will be moved 
		 * Order of atoms for improper as in the CHARMM IC:
		 *
		 *    1  4      Rotate 4 with respect to 1 on the
		 *     \/       2-3 axis.  
		 *  2--3
		 *
		 *******************************************************/
		bool setBondDistance(Atom & _atom1, Atom & _atom2, double _distance);
		bool setBondAngle(Atom & _atom1, Atom & _atom2, Atom & _atom3, double _angleDegrees);
		bool setDihedral(Atom & _atom1, Atom & _atom2, Atom & _atom3, Atom & _atom4, double _angleDegrees, bool _strict=false);
		bool setImproper(Atom & _atom1, Atom & _atom2, Atom & _atom3, Atom & _atom4, double _angleDegrees);



		Matrix createBasisTransformation(std::vector<std::vector<double> > &_basis1, std::vector<std::vector<double> > &_basis2);


		SphericalPoint transform(CartesianPoint &_p1);


		/*****************************************************
		 *  The object can store the history of all transofmations
		 *  as applied to a local frame.  If the history is saved,
		 *  the result of all recorded conformations can be applied
		 *  at once to any other atom or atom std::vector
		 *****************************************************/
		void setStoreTransformHistory(bool _flag);
		bool getStoreTransformHistory() const;
		void resetHistory();
		void applyHistory(Atom & _atom);
		void applyHistory(AtomPointerVector & _atoms);

		Matrix getLastRotationMatrix() const;
		CartesianPoint getLastTranslation() const;

		/*
		/ *******************************************************
                 * Functions to allow grid search (from Cinque Soto)
                 * RotatePdbAboutZYX rotates an AtomPointerVector and the three
                 * localAxes (centered on the origin) based on the input
                 * degree rotations in local Z, Y, and X.  Center is
                 * not changed.
                 * TranslateRigidBody translates AtomPointerVector and center
                 * my a specified distance.
                 ******************************************************* /
		void RotatePdbAboutZYX(AtomPointerVector & _theAtoms, CartesianPoint & _center, CartesianPoint & _localZ, CartesianPoint & _localY, CartesianPoint & _localX, double _RotationAlongZ, double _RotationAlongY, double _RotationAlongX);
		void TranslateRigidBodyPdbResidue(AtomPointerVector & _theAtoms, CartesianPoint & _center, CartesianPoint & _TranslationVector, double _TranslationAmount);
		*/

		// Move only the current or all coors?
		// if false only the current coors are moved
		void setTransformAllCoors(bool _flag); 
		bool getTransformAllCoors() const; 

	private:

	//	void findLinkedAtoms(Atom * _pAtom, const std::map<Atom*, bool> & _excluded, std::map<Atom*, bool> & _list);

		void translateAtom(Atom & _atom, const CartesianPoint & _p);
		void rotateAtom(Atom & _atom, const Matrix & _rotMatrix, const CartesianPoint & _rotCenter=CartesianPoint(0.0, 0.0, 0.0));
		bool alignAtom(Atom & _atom, const CartesianPoint & _target, const CartesianPoint & _rotCenter);
		bool orientAtom(Atom & _atom, const CartesianPoint & _target, const CartesianPoint & _axis1, const CartesianPoint & _axis2);

		Quaternion q;
		
		Matrix lastRotMatrix;
		CartesianPoint lastTranslation;
		bool saveHistory_flag;

		std::map<std::string, CartesianPoint> frame;

		bool transformAllCoors_flag; // if false only the current coors are moved
};

// INLINE FUNCTIONS
inline void Transforms::setStoreTransformHistory(bool _flag) {saveHistory_flag = _flag;}
inline bool Transforms::getStoreTransformHistory() const {return saveHistory_flag;}
inline void Transforms::resetHistory() {frame["O"] = CartesianPoint(0.0, 0.0, 0.0); frame["X"] = CartesianPoint(1.0, 0.0, 0.0); frame["Y"] = CartesianPoint(0.0, 1.0, 0.0);}
inline Matrix Transforms::getLastRotationMatrix() const {return lastRotMatrix;}
inline CartesianPoint Transforms::getLastTranslation() const {return lastTranslation;}
inline void Transforms::setTransformAllCoors(bool _flag) {transformAllCoors_flag = _flag;}
inline bool Transforms::getTransformAllCoors() const {return transformAllCoors_flag;} 
}

#endif
