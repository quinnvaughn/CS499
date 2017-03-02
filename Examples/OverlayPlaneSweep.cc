extern "C"{
  void overlayPlaneSweep( const halfsegment r1[], int r1Size,
												const halfsegment r2[], int r2Size,
												vector<halfsegment>& result )
{
	halfsegment currSeg, maxSeg, tmpSeg, *tmpSegPtr;
	maxSeg.dx = maxSeg.dy = maxSeg.sx = maxSeg.sy = std::numeric_limits<double>::max();
	double eventX, eventY;
	avl_table* activeList = avl_create( avlHsegActiveListCompare, &eventX, NULL );// active list avl tree
	avl_table* discoveredSegs = avl_create( avlHsegCompare, NULL, NULL ); // discovered segs
	avl_traverser discoveredTraverser;
	avl_traverser alTraverser, alAboveTraverser, alBelowTraverser;
	halfsegment* alHsegPtr = NULL, *alInsertPtr = NULL, *alHsegAbovePtr = NULL, *alHsegBelowPtr = NULL;
	vector< halfsegment > brokenSegs;
	bool colinearIntersection;
	avl_t_init( &discoveredTraverser, discoveredSegs );
	avl_t_init( &alTraverser, activeList );
	avl_t_init( &alAboveTraverser, activeList );
	avl_t_init( &alBelowTraverser, activeList );

	int r1Pos = 0;
	int r2Pos = 0;
	int segSource;
	while( discoveredSegs->avl_count > 0 || r1Pos < r1Size || r2Pos < r2Size ) {
		// get the next seg
		// next seg is the least seg from r1, r2, and the discoveredSeg tree (event queue)
		currSeg = maxSeg;
		if( r1Pos < r1Size ) {
			currSeg = r1[r1Pos];
			segSource = 1;
		}
		if( r2Pos < r2Size &&  r2[r2Pos]< currSeg ) {
			currSeg = r2[r2Pos];
			segSource = 2;
		}
		tmpSegPtr = static_cast<halfsegment*>(avl_t_first( &discoveredTraverser, discoveredSegs ) );
		if( tmpSegPtr != NULL ) {
			tmpSeg = *tmpSegPtr;
			if( tmpSeg < currSeg || tmpSeg == currSeg  ){
				currSeg = tmpSeg;
				segSource = 3;
			}
		}
		// remove the next seg from its source
		if( segSource == 3 ) {
			avl_delete( discoveredSegs, &currSeg );
			delete tmpSegPtr;
			tmpSegPtr = NULL;
		}
		else if( segSource == 2 ){
			r2Pos++;
		}
		else {
			r1Pos++;
		}

		// set current event point.
		// the activeList compare function uses eventX as its |param| argument
		eventX = currSeg.dx;
		eventY = currSeg.dy;

		// If curr is a left seg, insert it and check for intersections with neighbors
		// Else it is a right seg, remove it and check its neighbors for intersections
		if( currSeg.isLeft( ) ) {
			// initialize the overlap labels
			currSeg.ola = currSeg.olb = -1;
			// insert the left seg
			// use avl_t_insert so we can get its neighbors.
			alInsertPtr = new halfsegment( currSeg );
			alHsegPtr = static_cast<halfsegment*>( avl_t_insert( &alTraverser, activeList, alInsertPtr ));
			// if a duplicate is in the active list, we get a pointer to the duplicate. So we need to update labels
			// if insert is successful, we get a pointer to the inserted item
			if( alHsegPtr != alInsertPtr ) {
				// We found a duplicate in active list.  update labels
				delete alInsertPtr;
				alInsertPtr = NULL;
				// NOTE: overlapping segs are a special case for labels.
				// NOTE: overlap labels are altered in the |breakHsegs()| as well,
				//       but that function is not called here (no need to break up segs
				//       if they are equal)
				alHsegPtr->ola = currSeg.la;
				alHsegPtr->olb = currSeg.lb;
			}
			else {
				bool needToRemoveCurr = false;
				// inserted successfully.  Need to get neighbors
				avl_t_copy( &alAboveTraverser, &alTraverser );
				avl_t_copy( &alBelowTraverser, &alTraverser );
				alHsegBelowPtr = static_cast<halfsegment*>( avl_t_prev( &alBelowTraverser));
				alHsegAbovePtr = static_cast<halfsegment*>( avl_t_next( &alAboveTraverser));

				// do intersections with above and below.  Update currSeg along the way
				brokenSegs.clear();

				// We have to deal with the below seg first because currSegs labels get changed
				// based on the below seg.  Once we begin dealing with the above seg, the labels
				// for the currSeg are already computed and will carry over into those calaculations
				if( alHsegBelowPtr != NULL ) {
					halfsegment belowSegCopy(*alHsegBelowPtr);

					// update labels:
					// NOTE: overlapping segs are a special case for labels.
					// NOTE: overlap labels are altered in the |breakHsegs()| as well,
					if( currSeg.regionID != belowSegCopy.regionID ) {
						// if below seg is from opposing region, set overlap labels
						if( belowSegCopy.dx != belowSegCopy.sx ) { // if seg is not vertical, use la (label above)
							currSeg.ola = currSeg.olb = belowSegCopy.la;
						} else { // if below seg is vertical, use lb (the label to the right)
							currSeg.ola = currSeg.olb = belowSegCopy.lb;
						}
					}
					else if(currSeg.regionID == belowSegCopy.regionID ) { // if below seg is from same region, just extend overlap labels
						currSeg.ola = currSeg.olb = belowSegCopy.ola;
						// commented code is for checking against vertical seg below from same regions
						// this should never happen for well formed regions based on hseg order
						// if( belowSegCopy.dx != belowSegCopy.sx ) { // if seg is not vertical, use ola (label above)
						//		currSeg.ola = currSeg.olb = belowSegCopy.ola;
						//  } else { // if below seg is vertical, use lb (the label to the right)
						//	currSeg.ola = currSeg.olb = belowSegCopy.olb;
						//	}
					}

					// Labels are now computed
					// Compute the segment intersections:
					if(breakHsegs(  *alHsegBelowPtr, currSeg, brokenSegs, colinearIntersection, false ) ){
						needToRemoveCurr = true;
						// remove below seg
						alHsegBelowPtr = static_cast<halfsegment*>( avl_delete( activeList, alHsegBelowPtr ) );
						delete alHsegBelowPtr;
						alHsegBelowPtr = NULL;
					}

				}
				// compute intersections with above seg:
				if( alHsegAbovePtr != NULL ) {
					if( breakHsegs(  *alHsegAbovePtr, currSeg, brokenSegs, colinearIntersection, false ) ) {
						needToRemoveCurr = true;
						// remove above seg
						alHsegAbovePtr = static_cast<halfsegment*>( avl_delete( activeList, alHsegAbovePtr ) );
						delete alHsegAbovePtr;
						alHsegAbovePtr = NULL;
					}
				}

				// Update the seg inserted into the active list this round
				// currSeg is the result of intersecting that seg with its neighbors
				*alHsegPtr = currSeg;

				// Insert all the broken up segs into thier various data structures
				insertBrokenSegsToActiveListAndDiscoveredQueue( brokenSegs,result,
																													discoveredSegs, activeList,
																													eventX, eventY );

			}
		}
		else {
			// This is a right halfsegment.
			// find its brother (left halfsegment) in the active list,
			//      remove it, and check its neighbors for intersections.
			currSeg = currSeg.getBrother();
			alHsegPtr = static_cast<halfsegment*>( avl_t_find( &alTraverser, activeList, &currSeg ));
			if( alHsegPtr != NULL ) {
				// We found the halfsegment in the active list.
				// Its possible we don't find one, since a right halfsegment may be in r1 or r2
				//   whose brother (left halfsegment) was broken due to an intersection
				result.push_back( *alHsegPtr );
				result.push_back( alHsegPtr->getBrother( ) );
				// The copy of the seg in the active list has the appropriate labels, so copy it over
				currSeg = *alHsegPtr;
				// find the neighbors
				avl_t_copy( &alAboveTraverser, &alTraverser );
				avl_t_copy( &alBelowTraverser, &alTraverser );
				alHsegBelowPtr = static_cast<halfsegment*>( avl_t_prev( &alBelowTraverser));
				alHsegAbovePtr = static_cast<halfsegment*>( avl_t_next( &alAboveTraverser));
				// delete the segment
				alHsegPtr = static_cast<halfsegment*>(avl_delete( activeList, alHsegPtr ));
				delete alHsegPtr;
				alHsegPtr = NULL;
				// if both neighbors are not null, then there is an above and below neighbor,
				//    check them for an intersection
				if( alHsegBelowPtr != NULL && alHsegAbovePtr != NULL ) {
					tmpSeg = *alHsegAbovePtr;
					brokenSegs.clear();
					if(breakHsegs(  *alHsegBelowPtr, tmpSeg, brokenSegs, colinearIntersection, true ) ){
						// remove below seg and above seg
						alHsegBelowPtr = static_cast<halfsegment*>( avl_delete( activeList, alHsegBelowPtr ) );
						delete alHsegBelowPtr;
						alHsegBelowPtr = NULL;
						alHsegAbovePtr = static_cast<halfsegment*>( avl_delete( activeList, alHsegAbovePtr ) );
						delete alHsegAbovePtr;
						alHsegAbovePtr = NULL;
						// add broken segs to output and discovered list
						insertBrokenSegsToActiveListAndDiscoveredQueue( brokenSegs,result,
																														discoveredSegs, activeList,
																														eventX, eventY );
					}
				}

			}
		}

	}
		// sort the result

		std::sort( result.begin(), result.end() );
	//	cerr << "ps result: " <<endl;
	//for( int i = 0; i <result.size(); i++ )
	//	if( result[i].isLeft() )
	//		cerr << result[i]<<endl;

}
}
